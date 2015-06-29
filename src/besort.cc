#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

#include "tclap/CmdLine.h"

//-------------------------------------------------------------------------
//
//
// Options
//
//
//-------------------------------------------------------------------------

using Number = unsigned int;

struct Options {
    
    Options(std::vector<std::string>& args);
    
    TCLAP::CmdLine cmd_line { "besort", ' ', "0.1rc1", true };
    
    TCLAP::ValueArg<Number> record_size {
        "s",                    // flag
        "record-size",          // name
        "record size in bytes", // description
        true,                   // required
        0,                     // value
        "number"                // type description
    };

    TCLAP::ValueArg<Number> key_size {
        "k",             // flag
        "key-size",      // name
        "key size in bytes",  // description
        true,            // required
        0,              // value
        "number"         // type description
    };

    TCLAP::ValueArg<Number> key_offset {
        "f",             // flag
        "key-offset",       // name
        "key offset inside a record (in bytes)",     // description
        true,            // required
        0,              // value
        "number"         // type description
    };

    TCLAP::ValueArg<Number> header_offset {
        "h",             // flag
        "header-offset",       // name
        "header offset in bytes",     // description
        false,           // required
        0,               // value
        "filename"       // type description
    };

    // -s or --schema
    TCLAP::ValueArg<std::string> buffer_size {
        "b",                     // flag
        "buffer-size",           // name
        "buffer size (no suffix or b=bytes, k=kilo-bytes, m=mega-bytes, g=giga-bytes",  // description
        false,                   // required
        "128m",                  // 128 megabytes is the default
        "size"                   // type description
    };

    TCLAP::ValueArg<std::string> input {
        "i",                     // flag
        "input-filename",        // name
        "input file name",       // description
        true,                    // required
        "",                       // value
        "string"                 // type description
    };

    TCLAP::ValueArg<std::string> output {
        "o",                     // flag
        "output-filename",        // name
        "output file name",       // description
        true,                    // required
        "",                       // value
        "string"                 // type description
    };

};

//------------------------------------------------------------------------------
// Options Impl.
//------------------------------------------------------------------------------

Options::Options(std::vector<std::string>& args)
{
    cmd_line.add(record_size);
    cmd_line.add(key_size);
    cmd_line.add(key_offset);
    cmd_line.add(header_offset);
    cmd_line.add(buffer_size);
    cmd_line.add(input);
    cmd_line.add(output);
    cmd_line.parse(args);
}



//------------------------------------------------------------------------------
// Forward Declarations
//------------------------------------------------------------------------------

struct BESort;
struct Chunk;

using  RecordCount  = std::uint64_t;
using  RecordSize   = std::uint64_t;
using  RecordOffset = std::uint64_t;
using  ChunkCount   = std::uint64_t;
using  FileOffset   = std::uint64_t;
using  FileSize     = std::uint64_t;
using  BufferSize   = std::uint64_t;
using  Index        = std::uint64_t;

// global variables

static RecordSize   RECORD_SIZE = 0; // statically set this variable
static RecordSize   KEY_SIZE    = 0;
static RecordOffset KEY_OFFSET  = 0;

//------------------------------------------------------------------------------
// Chunk
//------------------------------------------------------------------------------

struct Chunk {
public:
    Chunk(BESort *parent, Index index, FileOffset begin, FileOffset end, RecordCount num_records);
    std::string filename() const;
public:
    BESort*     parent { nullptr };
    Index       index { 0 };
    
    FileOffset  begin { 0 };
    FileOffset  end   { 0 };
    RecordCount num_records { 0 };
};

//------------------------------------------------------------------------------
// Chunk Impl.
//------------------------------------------------------------------------------

Chunk::Chunk(BESort *parent,
             Index   index,
             FileOffset begin,
             FileOffset end,
             RecordCount num_records):
parent(parent),
index(index),
begin(begin),
end(end),
num_records(num_records)
{}


std::string Chunk::filename() const {
    std::stringstream ss;
    ss << "/tmp/__besort" << index;
    return ss.str();
}

//------------------------------------------------------------------------------
// BESort
//------------------------------------------------------------------------------

/*! \brief Binary External Sort
 */
struct BESort {
public:    
    BESort(const std::string& input_filename,
            const std::string& output_filename,
            FileOffset         header_offset,
            RecordSize         record_size,
            RecordSize         key_size,
            RecordOffset       key_offset,
            BufferSize         buffer_size);
public:    

    struct {
        std::string   input_filename;
        std::string   output_filename;
        FileOffset    header_offset;
        RecordSize    record_size;
        RecordSize    key_size;
        RecordOffset  key_offset;
        BufferSize    buffer_size;
    } input;

    struct {
        FileSize file_size;
        RecordCount num_records;
        RecordCount records_per_chunk;
        ChunkCount  num_chunks;
    } computed;
    
    std::vector<std::unique_ptr<Chunk>> chunks;
};

//------------------------------------------------------------------------------
// BESort Impl.
//------------------------------------------------------------------------------

BESort::BESort(
               const std::string& input_filename,
               const std::string& output_filename,
               FileOffset   header_offset,
               RecordSize   record_size,
               RecordSize   key_size,
               RecordOffset key_offset,
               BufferSize buffer_size)
{
    // input fields
    input.input_filename = input_filename;
    input.output_filename = output_filename;
    input.header_offset = header_offset;
    input.record_size = record_size;
    input.key_size    = key_size;
    input.key_offset = key_offset;
    input.buffer_size = buffer_size;
    
    // computed fields
    
    // go to header
    std::ifstream is(input.input_filename);
    if (!is)
        throw std::runtime_error("couldn't open file");

    // file size
    is.seekg(0, is.end);
    computed.file_size = is.tellg();
    
    // number of records (integer division here)
    computed.num_records = (computed.file_size - input.header_offset)/input.record_size;
    
    // chunks will be of size
    computed.records_per_chunk = input.buffer_size / (std::uint64_t) input.record_size;

    // number of chunks
    computed.num_chunks = computed.num_records / computed.records_per_chunk;

    // last chunk might need
    computed.num_chunks += (computed.num_records & computed.records_per_chunk) ? 1 : 0;
    

    // prepare chunks
    auto remaining = computed.num_records;
    auto offset    = input.header_offset;

    chunks.reserve(computed.num_chunks);
    
    
    // set global variables
    KEY_OFFSET  = input.key_offset;
    KEY_SIZE    = input.key_size;
    RECORD_SIZE = input.record_size;
    
    
    auto index = 0;
    while (remaining > 0) {
        auto r = (remaining < computed.records_per_chunk) ? remaining : computed.records_per_chunk;
        auto s = (r * input.record_size);
        auto e = offset + s;
        chunks.push_back(std::unique_ptr<Chunk>(new Chunk(this, index++, offset, e, r)));
        offset = e;
        remaining -= r;
    }
    
}






//------------------------------------------------------------------------------
// Record Iterator
//------------------------------------------------------------------------------

using RecordDiffType = int64_t;

struct Record {
    Record() = default;
    
    Record(const Record& other);
    Record& operator=(const Record &other);
    
    bool operator<(const Record& other) const;
    bool operator==(const Record& other) const;
};

Record::Record(const Record& other) {
    auto pa = (char*) this;
    auto pb = (const char*) &other;
    auto eb = pb + RECORD_SIZE;
    std::copy(pb,eb,pa);
}

Record& Record::operator=(const Record &other) {
    auto pa = (char*) this;
    auto pb = (const char*) &other;
    auto eb = pb + RECORD_SIZE;
    std::copy(pb,eb,pa);
    return *this;
}

bool Record::operator<(const Record& other) const {
    if (KEY_SIZE <= 8) {
        uint64_t a = 0;
        uint64_t b = 0;
        auto pa = (char*) this   + KEY_OFFSET;
        auto pb = (char*) &other + KEY_OFFSET;
        std::copy(pa, pa + KEY_SIZE, (char*) &a);
        std::copy(pb, pb + KEY_SIZE, (char*) &b);
        return a < b;
    }
    else {
        throw std::runtime_error("not implemented yet");
    }
}

bool Record::operator==(const Record& other) const {
    if (KEY_SIZE == 8) {
        uint64_t a = 0;
        uint64_t b = 0;
        auto pa = (char*) this   + KEY_OFFSET;
        auto pb = (char*) &other + KEY_OFFSET;
        std::copy(pa, pa + KEY_SIZE, (char*) &a);
        std::copy(pb, pb + KEY_SIZE, (char*) &b);
        return a < b;
        return a == b;
    }
    else {
        throw std::runtime_error("not implemented yet");
    }
}

static inline Record* rp_add(const Record* p, RecordDiffType n) {
    return (Record*) ((char*) p + n * RECORD_SIZE);
}

//------------------------------------------------------------------------------
// Record Iterator
//------------------------------------------------------------------------------

struct RecordIterator:
public std::iterator<std::random_access_iterator_tag, Record>
{
public:
    using type              = std::iterator<std::random_access_iterator_tag, Record>;
    using iterator_category = std::random_access_iterator_tag;
    using value_type        = type::value_type;
    using difference_type   = type::difference_type;
    using reference         = type::reference;
    using pointer           = type::pointer;
    
    RecordIterator() {};
    
    inline RecordIterator(const RecordIterator& r) : ptr(&(*r))
    {}
    
    inline RecordIterator(pointer p) : ptr(p) {}
    
    inline RecordIterator& operator=(const RecordIterator& r)
    { ptr = &(*r); return *this; }
    
    inline RecordIterator& operator++()  // PREFIX
    { ptr = rp_add(ptr,1); return *this; }

    inline RecordIterator& operator--()  // PREFIX
    { ptr = rp_add(ptr,-1); return *this; }

    inline RecordIterator operator++(int)  // POSTFIX
    {
        auto result = *this;
        ptr = rp_add(ptr,1);
        return result;
    }

    RecordIterator operator--(int)  // POSTFIX
    {
        auto result = *this;
        ptr = rp_add(ptr,-1);
        return result;
    }
    
    RecordIterator operator+(const difference_type& n) const
    { return RecordIterator(rp_add(ptr,n)); }
    
    RecordIterator& operator+=(const difference_type& n)
    { ptr = (Record*) ((char*) ptr + n * RECORD_SIZE); return *this; }
    
    RecordIterator operator-(const difference_type& n) const
    { return RecordIterator(pointer(rp_add(ptr,-n))); }
    
    RecordIterator& operator-=(const difference_type& n) {
        ptr = rp_add(ptr,-n); return *this;
    }

    bool operator==(const RecordIterator& other) const {
        return ptr == other.ptr;;
    }

    reference operator*() const { return *ptr; }
    
    pointer operator->() const { return ptr; }
    
    reference operator[](const difference_type& n) const { return *(rp_add(ptr,n)); }

protected:
    Record* ptr;
};


//------------------------------------------------------------------------------
// LoadChunkSortAndSave
//------------------------------------------------------------------------------

struct LoadChunkSortAndSave {
public:
    LoadChunkSortAndSave(Chunk &chunk);
    void run();
public:
    Chunk *chunk_p { nullptr };
    std::vector<char> buffer;
};


LoadChunkSortAndSave::LoadChunkSortAndSave(Chunk &chunk):
    chunk_p(&chunk)
{}

void LoadChunkSortAndSave::run() {
    
    //
    auto &chunk = *chunk_p;
    buffer.resize(chunk.end - chunk.begin); //
    
    
    { // read buffer
        std::ifstream is(chunk.parent->input.input_filename);
        is.seekg(chunk.begin);
        is.read(&buffer[0], buffer.size());
    }
    
    
    auto p     = &buffer[0];
    auto begin = (Record*) p;
    auto end   = (Record*) (p + buffer.size()) ;
    
    
    // in memory sort
    std::sort(begin, end);
    
    //
    std::ofstream of(chunk.filename());
    of.write(&buffer[0],buffer.size());
}

//------------------------------------------------------------------------------
// Queue
//------------------------------------------------------------------------------

struct Queue {
public:

    Queue(Chunk* chunk, RecordCount buffer_record_capacity);

    const Record* front() const;
    
    Record* front();

    bool operator<(const Queue& other) const;

    Record* next();
    
public:
    std::vector<char> buffer;
    Chunk*            chunk;
    RecordCount       buffer_record_capacity;
    RecordCount       remaining;
    std::ifstream     is;
    RecordIterator    it;
    RecordIterator    it_end;
};

//------------------------------------------------------------------------------
// Queue Impl.
//------------------------------------------------------------------------------

Queue::Queue(Chunk* chunk, BufferSize buffer_record_capacity):
    chunk(chunk),
    buffer_record_capacity(buffer_record_capacity),
    is(chunk->filename())
{
    remaining = chunk->num_records;
    
    auto batch = remaining > buffer_record_capacity ? buffer_record_capacity : remaining;
    buffer.resize(batch);
    
    auto p = &buffer[0];
    is.read(p, batch * chunk->parent->input.record_size);
    
    it     = RecordIterator((Record*) &buffer[0]);
    it_end = it + batch;
}

const Record* Queue::front() const{
    if (remaining == 0)
        return nullptr;
    else
        return &(*it);
}

Record* Queue::front() {
    if (remaining == 0)
        return nullptr;
    else
        return &(*it);
}

bool Queue::operator<(const Queue& other) const{
    return front() > other.front();
}

Record* Queue::next() {
    
    if (remaining == 0)
        return nullptr;
        
    auto result = front();
    
    ++it;
    --remaining;
    
    // advance
    if (it == it_end && remaining > 0) {
        auto batch = remaining > buffer_record_capacity ? buffer_record_capacity : remaining;
        buffer.resize(batch);
        
        auto p = &buffer[0];
        is.read(p, batch * chunk->parent->input.record_size);
        
        it     = RecordIterator((Record*) &buffer[0]);
        it_end = it + batch;
    }
    
    return result;
}


//------------------------------------------------------------------------------
// MultiWayMerge
//------------------------------------------------------------------------------

struct MultiWayMerge {
public:
    MultiWayMerge(BESort *parent, Index index, const std::vector<Chunk*> &chunks, RecordCount buffer_record_capacity);
    void run();
    std::string filename() const;
public:
    BESort *parent { nullptr };
    Index  index;
    std::vector<std::unique_ptr<Queue>> queues;
    RecordCount buffer_record_capacity;
};


std::string MultiWayMerge::filename() const {
    std::stringstream ss;
    ss << "__besort_mway_" << index;
    return ss.str();
}


MultiWayMerge::MultiWayMerge(BESort *parent, Index index, const std::vector<Chunk*> &chunks, RecordCount buffer_record_capacity):
parent(parent),
index(index),
buffer_record_capacity(buffer_record_capacity)
{
    
    // initialize queues: one for each chunk
    for (auto c: chunks) {
        queues.push_back(std::unique_ptr<Queue>{ new Queue(c, buffer_record_capacity) });
    }
    
    std::make_heap(queues.begin(),queues.end(),[](const std::unique_ptr<Queue>& a, const std::unique_ptr<Queue>& b) {
        auto &qa = *a.get();
        auto &qb = *b.get();
        return qa < qb;
    });
}

void MultiWayMerge::run() {

    std::ofstream os(filename());
    while (queues.size()) {
        auto record = queues.front()->next();
        os.write((char*) record, parent->input.record_size);
        
        // empty
        auto next_front = queues.front()->front();
        if (!next_front) {
            std::pop_heap(queues.begin(),queues.end(),[](const std::unique_ptr<Queue>& a, const std::unique_ptr<Queue>& b) {
                auto &qa = *a.get();
                auto &qb = *b.get();
                return qa < qb;
            });
            queues.pop_back();
        }
        else if (queues.size() > 1) {
            bool swap_queues = false;
            // children
//
            auto &q_root  = *queues.begin()->get();
//            auto &q_child1 = *queues[1].get();
//            if (q_child1 < q_root) {
//                swap_queues = true;
//            }
//            else if (queues.size() > 2) {
//                auto &q_child2 = *queues[2].get();
//                if (q_child2 < q_root) {
//                    swap_queues = true;
//                }
//            }
            if (swap_queues) {
                std::pop_heap(queues.begin(),queues.end(),[](const std::unique_ptr<Queue>& a, const std::unique_ptr<Queue>& b) {
                    auto &qa = *a.get();
                    auto &qb = *b.get();
                    return qa < qb;
                });
                std::push_heap(queues.begin(),queues.end(),[](const std::unique_ptr<Queue>& a, const std::unique_ptr<Queue>& b) {
                    auto &qa = *a.get();
                    auto &qb = *b.get();
                    return qa < qb;
                });
            }
        }
        
        
        
    }

}



//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------

int main(int argc, char** argv) {

    // load options
    std::vector<std::string> params(argv, argv+argc);
    Options options(params);

    // print program name
    std::cerr << options.cmd_line.getMessage() << "  v." << options.cmd_line.getVersion() << std::endl;

    // figure out file size
    std::ifstream is(options.input.getValue());

    if (!is)
        throw std::runtime_error("Could not open input file");


    auto parse_size = [](const std::string& st) -> std::uint64_t {
        if (!st.size())
            return 0;

        uint64_t multiplier = 0;
        
        char unit = std::tolower(st.back());
        
        if (unit == 'b') {
            multiplier = 1;
        }
        else if (unit == 'k') {
            multiplier = 1 << 10;
        }
        else if (unit == 'm') {
            multiplier = 1 << 20;
        }
        else if (unit == 'g') {
            multiplier = 1 << 30;
        }
        
        auto num = std::stoull(std::string(st.begin(),
                                           st.end() + ((multiplier > 0) ? -1 : 0) ) );
        
        if (multiplier == 0)
            multiplier = 1;
        
        return multiplier * num;
    };
    

    //
    // setup BESort object
    //
    
    BESort besort {
        options.input.getValue(),
        options.output.getValue(),
        options.header_offset.getValue(),
        options.record_size.getValue(),
        options.key_size.getValue(),
        options.key_offset.getValue(),
        parse_size(options.buffer_size.getValue())
    };
    
    //
    //
    //
    for (auto &it: besort.chunks) {
        LoadChunkSortAndSave task(*it.get());
        task.run();
    }
    
    
}




    
