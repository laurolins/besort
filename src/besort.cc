#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstdio>

#include "tclap/CmdLine.h"

#include "thread_pool.hh"

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

    TCLAP::ValueArg<int> threads {
        "t",                     // flag
        "threads",               // name
        "output file name",      // description
        false,                   // required
        4,                       // value
        "number"                 // type description
    };

    TCLAP::ValueArg<int> multiway_merge {
        "m",                     // flag
        "mulit-way-merge",               // name
        "multi way merge parameter (default == 2)",      // description
        false,                   // required
        2,                       // value
        "number"                 // type description
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
    cmd_line.add(threads);
    cmd_line.add(multiway_merge);
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
using  ChunkLevel   = std::uint64_t; // level 0 are leaves, level 1 are chunks composed of level 0 chunks
using  ChunkCount   = std::uint64_t; // level 0 are leaves, level 1 are chunks composed of level 0 chunks

// global variables

static RecordSize   RECORD_SIZE = 0; // statically set this variable
static RecordSize   KEY_SIZE    = 0;
static RecordOffset KEY_OFFSET  = 0;

//------------------------------------------------------------------------------
// Chunk
//------------------------------------------------------------------------------

struct Chunk {
public:
    Chunk(const std::vector<Chunk*>& children);
    Chunk(BESort *besort, Chunk* parent, Index index, FileOffset begin, FileOffset end, RecordCount num_records, ChunkLevel level);
    
    /*!
     * Returns true if this non-leaf chunk is ready to be merged
     * Updates num sorted children and when it reaches number of children
     * returns true
     */
    bool aChildWasSorted();
    
    std::string filename() const;
public:
    BESort*     besort { nullptr };
    Chunk*      parent { nullptr };
    Index       index { 0 };
    
    ChunkLevel  level { 0 };
    FileOffset  begin { 0 };
    FileOffset  end   { 0 };
    RecordCount num_records         { 0 };

    std::vector<Chunk*> children;
    ChunkCount          num_sorted_children { 0 };
};

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
            BufferSize         buffer_size,
            int                num_threads,
            int                multiway_merge);
private:
    
    void _run();
    
public:
    
    
    
public:

    struct {
        std::string   input_filename;
        std::string   output_filename;
        FileOffset    header_offset;
        RecordSize    record_size;
        RecordSize    key_size;
        RecordOffset  key_offset;
        BufferSize    buffer_size;
        int           num_threads;
        int           multiway_merge;
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
// Record Iterator
//------------------------------------------------------------------------------

using RecordDiffType = int64_t;

struct Record {
    Record() = default;
    
    Record(const Record& other);
    Record& operator=(const Record &other);

    Record(Record&& other);
    Record& operator=(Record &&other);
    
    uint64_t key() const; // if KEY_SIZE <= 8 bytes key() works
    
    bool operator<(const Record& other) const;
    bool operator>(const Record& other) const;
    bool operator==(const Record& other) const;
};

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
    
    RecordIterator() = default;
    
    inline RecordIterator(const RecordIterator& r);
    inline RecordIterator(pointer p);
    inline RecordIterator& operator=(const RecordIterator& r);
    inline RecordIterator& operator++();  // PREFIX
    inline RecordIterator& operator--();  // PREFIX
    inline RecordIterator operator++(int);  // POSTFIX
    
    auto operator-(const RecordIterator& other) const -> difference_type;

    bool operator<(const RecordIterator& other) const;
    bool operator>=(const RecordIterator& other) const;
    bool operator>(const RecordIterator& other) const;
    bool operator<=(const RecordIterator& other) const;   
    
    RecordIterator operator--(int);  // POSTFIX
    RecordIterator operator+(const difference_type& n) const;
    RecordIterator& operator+=(const difference_type& n);
    
    RecordIterator operator-(const difference_type& n) const;
    RecordIterator& operator-=(const difference_type& n);

    bool operator==(const RecordIterator& other) const;
    bool operator!=(const RecordIterator& other) const;

    reference operator*() const;
    pointer operator->() const;
    reference operator[](const difference_type& n) const;

protected:
    Record* ptr { nullptr };
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
    std::vector<char> next_result;
    Chunk*            chunk;
    RecordCount       buffer_record_capacity;
    RecordCount       remaining;
    std::ifstream     is;
    RecordIterator    it;
    RecordIterator    it_end;
};


//------------------------------------------------------------------------------
// MultiWayMerge
//------------------------------------------------------------------------------

struct MultiWayMerge {
public:
    MultiWayMerge() = default;
    MultiWayMerge(Chunk *chunk);
    void run();
public:
    Chunk  *chunk;
    std::vector<std::unique_ptr<Queue>> queues;
};








































//------------------------------------------------------------------------------
// Chunk Impl.
//------------------------------------------------------------------------------

Chunk::Chunk(BESort* besort,
             Chunk*  parent,
             Index   index,
             FileOffset begin,
             FileOffset end,
             RecordCount num_records,
             ChunkLevel level):
besort(besort),
parent(parent),
index(index),
begin(begin),
end(end),
num_records(num_records),
level(level)
{}

Chunk::Chunk(const std::vector<Chunk*>& children):
    children(children)
{
    if (children.size() < 1)
        throw std::runtime_error("Higher level chunk should have at least two children");
    
    besort = children[0]->besort;
    level  = children[0]->level + 1;
    index  = children[0]->index;
    begin  = children[0]->begin;
    end    = children[0]->end;

    for (auto c: children) {
        index        = std::min(index, c->index);
        begin        = std::min(begin, c->begin);
        end          = std::max(end,   c->end);
        num_records += c->num_records;
        c->parent = this;
    }
    
    num_sorted_children = 0;
}

bool Chunk::aChildWasSorted() {
    ++num_sorted_children;
    return num_sorted_children == children.size();
}

std::string Chunk::filename() const {
    std::stringstream ss;
    ss << "/tmp/__besort_chunk_l" << level << "_i" << index << ".txt";
    return ss.str();
}


//------------------------------------------------------------------------------
// Record Impl.
//------------------------------------------------------------------------------

namespace std {
    template<>
    void swap(Record& ra, Record& rb)
    {
        auto pa = (char*) &ra;
        auto pb = (char*) &rb;
        for (int i=0;i<RECORD_SIZE;++i) {
            std::swap(*pa, *pb);
            ++pa;
            ++pb;
        }
    }
}

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

Record::Record(Record&& other) {
    auto pa = (char*) this;
    auto pb = (const char*) &other;
    auto eb = pb + RECORD_SIZE;
    std::copy(pb,eb,pa);
}

Record& Record::operator=(Record &&other) {
    auto pa = (char*) this;
    auto pb = (const char*) &other;
    auto eb = pb + RECORD_SIZE;
    std::copy(pb,eb,pa);
    return *this;
}

uint64_t Record::key() const {
    if (KEY_SIZE <= 8) {
        uint64_t a = 0;
        auto pa = (char*) this   + KEY_OFFSET;
        std::copy(pa, pa + KEY_SIZE, (char*) &a);
        return a;
    }
    else {
        throw std::runtime_error("not implemented yet");
    }
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

bool Record::operator>(const Record& other) const {
    if (KEY_SIZE <= 8) {
        uint64_t a = 0;
        uint64_t b = 0;
        auto pa = (char*) this   + KEY_OFFSET;
        auto pb = (char*) &other + KEY_OFFSET;
        std::copy(pa, pa + KEY_SIZE, (char*) &a);
        std::copy(pb, pb + KEY_SIZE, (char*) &b);
        return a > b;
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

inline RecordIterator::RecordIterator(const RecordIterator& r) : ptr(&(*r))
{}

inline RecordIterator::RecordIterator(pointer p) : ptr(p) {}

inline RecordIterator& RecordIterator::operator=(const RecordIterator& r)
{ ptr = &(*r); return *this; }

inline RecordIterator& RecordIterator::operator++()  // PREFIX
{ ptr = rp_add(ptr,1); return *this; }

inline RecordIterator& RecordIterator::operator--()  // PREFIX
{ ptr = rp_add(ptr,-1); return *this; }

inline RecordIterator RecordIterator::operator++(int)  // POSTFIX
{
    auto result = *this;
    ptr = rp_add(ptr,1);
    return result;
}

RecordIterator RecordIterator::operator--(int)  // POSTFIX
{
    auto result = *this;
    ptr = rp_add(ptr,-1);
    return result;
}

RecordIterator RecordIterator::operator+(const difference_type& n) const
{ return RecordIterator(rp_add(ptr,n)); }

RecordIterator& RecordIterator::operator+=(const difference_type& n)
{ ptr = (Record*) ((char*) ptr + n * RECORD_SIZE); return *this; }

RecordIterator RecordIterator::operator-(const difference_type& n) const
{ return RecordIterator(pointer(rp_add(ptr,-n))); }

auto RecordIterator::operator-(const RecordIterator& other) const -> difference_type
{ return ((difference_type)ptr - (difference_type)other.ptr)/RECORD_SIZE; }


RecordIterator& RecordIterator::operator-=(const difference_type& n) {
    ptr = rp_add(ptr,-n); return *this;
}

bool RecordIterator::operator==(const RecordIterator& other) const {
    return ptr == other.ptr;
}

bool RecordIterator::operator<(const RecordIterator& other) const {
    return ptr < other.ptr;
}

bool RecordIterator::operator>(const RecordIterator& other) const {
    return ptr > other.ptr;
}

bool RecordIterator::operator<=(const RecordIterator& other) const {
    return ptr <= other.ptr;
}

bool RecordIterator::operator>=(const RecordIterator& other) const {
    return ptr >= other.ptr;
}

bool RecordIterator::operator!=(const RecordIterator& other) const {
    return ptr != other.ptr;
}


auto RecordIterator::operator*() const -> reference
{ return *ptr; }

auto RecordIterator::operator->() const -> pointer
{ return ptr; }

auto RecordIterator::operator[](const difference_type& n) const -> reference {
    return *(rp_add(ptr,n));
}

//------------------------------------------------------------------------------
// LoadChunkSortAndSave Impl.
//------------------------------------------------------------------------------

LoadChunkSortAndSave::LoadChunkSortAndSave(Chunk &chunk):
chunk_p(&chunk)
{}

void LoadChunkSortAndSave::run() {
    
    //
    auto &chunk = *chunk_p;
    buffer.resize(chunk.end - chunk.begin); //
    
    
    { // read buffer
        std::ifstream is(chunk.besort->input.input_filename);
        is.seekg(chunk.begin);
        is.read(&buffer[0], buffer.size());
    }
    
    
    auto p     = &buffer[0];
    RecordIterator begin{(Record*) p};
    RecordIterator end  {(Record*) (p + buffer.size())};
    
    // std::cerr << std::distance(begin,end) << std::endl;
    
    // in memory sort
    std::sort(begin, end);
    
    //
    std::ofstream of(chunk.filename());
    of.write(&buffer[0],buffer.size());
}


//------------------------------------------------------------------------------
// Queue Impl.
//------------------------------------------------------------------------------

Queue::Queue(Chunk* chunk, BufferSize buffer_record_capacity):
    chunk(chunk),
    buffer_record_capacity(buffer_record_capacity),
    is(chunk->filename())
{
//    std::cerr << "Queue(chunk): " << chunk->filename() << std::endl;
    
    remaining = chunk->num_records;

//    std::cerr << "Queue(remaining): " << remaining << std::endl;

    auto batch = remaining > buffer_record_capacity ? buffer_record_capacity : remaining;
    auto batch_bytes = batch * chunk->besort->input.record_size;
    buffer.resize(batch_bytes);

    next_result.resize(chunk->besort->input.record_size);
    
    auto p = &buffer[0];
    is.read(p, batch_bytes);
    
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
    auto &a = *front();
    auto &b = *other.front();
    return a > b; // top priority is reversed on push_heap pop_heap etc
}

Record* Queue::next() {
    
    if (remaining == 0)
        return nullptr;
    
    Record* result = (Record*) &next_result[0];
    *result = *front();
    
    ++it;
    --remaining;
    
    // advance
    if (it == it_end && remaining > 0) {
        auto batch = remaining > buffer_record_capacity ? buffer_record_capacity : remaining;
        auto batch_bytes = batch * chunk->besort->input.record_size;
        buffer.resize(batch_bytes);
        
        auto p = &buffer[0];
        is.read(p, batch_bytes);
        
        it     = RecordIterator((Record*) &buffer[0]);
        it_end = it + batch;
    }
    
    return result;
}













//------------------------------------------------------------------------------
// MultiWayMerge Impl.
//------------------------------------------------------------------------------

MultiWayMerge::MultiWayMerge(Chunk* chunk):
chunk(chunk)
{

    // record capacity of each queue is given by
    auto &input = chunk->besort->input;
    auto buffer_record_capacity = input.buffer_size/(input.record_size * input.multiway_merge);
    if (buffer_record_capacity == 0)
        buffer_record_capacity = 1;
    

    // initialize queues: one for each chunk
    for (auto c: chunk->children) {
        queues.push_back(std::unique_ptr<Queue>{ new Queue(c, buffer_record_capacity) });
    }
    
    std::make_heap(queues.begin(),queues.end(),[](const std::unique_ptr<Queue>& a, const std::unique_ptr<Queue>& b) {
        auto &qa = *a.get();
        auto &qb = *b.get();
        return qa < qb;
    });
}

void MultiWayMerge::run() {
    
    
    auto compare_heaps = [](const std::unique_ptr<Queue>& a, const std::unique_ptr<Queue>& b) {
        auto &qa = *a.get();
        auto &qb = *b.get();
        return qa < qb;
    };

    std::ofstream os(chunk->filename());
    while (queues.size()) {
        auto record = queues.front()->next();
        os.write((char*) record, chunk->besort->input.record_size);
        
        // empty
        auto next_front = queues.front()->front();
        if (!next_front) {
            std::pop_heap(queues.begin(),queues.end(),compare_heaps);
            queues.pop_back();
        }
        else if (queues.size() > 1) {
            bool swap_queues = false;
            // children


            // std::cerr << "comparing q[0]=" << queues[0].get()->front()->key() << "  q[1]=" << queues[1].get()->front()->key() << std::endl;
            
            // auto &qroot = *queues.begin()->get();
            if (compare_heaps(queues[0], queues[1])) {
                swap_queues = true;
            }
            else if (queues.size() > 2 && compare_heaps(queues[0], queues[2])) {
                swap_queues = true;
            }

            
            if (swap_queues) {
                std::pop_heap(queues.begin(),queues.end(),compare_heaps);
                std::push_heap(queues.begin(),queues.end(),compare_heaps);
            }
        }
    }
}

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
               BufferSize buffer_size,
               int        num_threads,
               int        multiway_merge)
{
    // input fields
    input.input_filename  = input_filename;
    input.output_filename = output_filename;
    input.header_offset   = header_offset;
    input.record_size     = record_size;
    input.key_size        = key_size;
    input.key_offset      = key_offset;
    input.buffer_size     = buffer_size;
    input.num_threads     = num_threads;
    input.multiway_merge  = multiway_merge;
    
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
    
    
    
    //
    // Create leaf chunks
    //
    auto prev_size = chunks.size();
    
    auto index = 0;
    while (remaining > 0) {
        auto r = (remaining < computed.records_per_chunk) ? remaining : computed.records_per_chunk;
        auto s = (r * input.record_size);
        auto e = offset + s;
        chunks.push_back(std::unique_ptr<Chunk>(new Chunk(this, nullptr, index++, offset, e, r, 0)));
        offset = e;
        remaining -= r;
    }
    
    std::cerr << "chunks in level 0: " << chunks.size() << std::endl;
    
    //
    // prepare chunks at other levels
    //
    auto level = 1;
    while (chunks.size() > prev_size + 1) {
        // at least two chunks need to be merged in the next level
        auto i = prev_size;
        auto e = std::min(prev_size + input.multiway_merge,chunks.size());

        prev_size = chunks.size(); // setup fpr next chunk size
        
        while (i < e) {
            std::vector<Chunk*> items(e-i);
            std::transform(chunks.begin()+i,chunks.begin()+e,items.begin(),[](const std::unique_ptr<Chunk> &it) {
                return it.get();
            });
            chunks.push_back(std::unique_ptr<Chunk>(new Chunk(items)));
            i = e;
            e = std::min(i + input.multiway_merge,prev_size);
        }

        std::cerr << "chunks in level " << level++ << ": " << (chunks.size()-prev_size) << std::endl;
    }
    
    _run();

}


void BESort::_run() {
    
    //
    // pipeline
    //
    // let worker thread solve LoadChunkSortAndSave
    //
    // k-way sort.
    //
    
    
    thread_pool::ThreadPool pool(input.num_threads);
    
    auto chunks_to_sort = chunks.size();
    std::mutex update_mutex;
//    std::vector<Chunk*> sorted_chunks;
    
    auto &besort = *this;
    
    std::function<void(Chunk*)> sort;
    sort = [&besort, &sort, &pool, &update_mutex, &chunks_to_sort](Chunk* chunk) {
        if (chunk->level == 0) {
            {
                std::lock_guard<std::mutex> guard(update_mutex);
                std::cerr << "sorting leaf chunk: " << chunk->filename() << std::endl;
            }
            LoadChunkSortAndSave task(*chunk);
            task.run();
            {
                std::lock_guard<std::mutex> guard(update_mutex);
                --chunks_to_sort;
                std::cerr << "finished sorting: " << chunk->filename() << std::endl;
                auto parent_chunk = chunk->parent;
                if (parent_chunk) {
                    bool process_parent_chunk = chunk->parent->aChildWasSorted();
                    if (process_parent_chunk) {
                        pool.enqueue([parent_chunk, &sort]() { sort(parent_chunk); });
                    }
                }
            }
        }
        else {
            {
                std::lock_guard<std::mutex> guard(update_mutex);
                std::cerr << "merging sorted chunks: " << chunk->filename() << std::endl;
            }

            MultiWayMerge multiway_merge { chunk };
            multiway_merge.run();
            
            // delete tmp files of children chunks
            for (auto c: chunk->children) {
                std::remove(c->filename().c_str());
            }

            {
                std::lock_guard<std::mutex> guard(update_mutex);
                --chunks_to_sort;
                std::cerr << "finished merging: " << chunk->filename() << std::endl;

                auto parent_chunk = chunk->parent;
                if (parent_chunk) {
                    bool process_parent_chunk = parent_chunk->aChildWasSorted();
                    if (process_parent_chunk) {
                        pool.enqueue([parent_chunk, &sort]{ sort(parent_chunk); });
                    }
                }
            }
        }
    };
    
    for (auto &it: chunks) {
        auto chunk = it.get();
        if (chunk->level == 0) {
            pool.enqueue([chunk, &sort] { sort(chunk); });
        }
    }
    
    while (chunks_to_sort > 0) {
        std::this_thread::sleep_for(std::chrono::milliseconds{200});
    }
    
    
    // rename file to output filename
    if (input.header_offset == 0){
        std::rename(chunks.back()->filename().c_str(),input.output_filename.c_str());
    }
    else {
        // copy header
        std::ofstream os(input.output_filename);
        
        { // copy header
            std::ifstream is(input.input_filename);
            std::size_t remaining = input.header_offset;
            std::vector<char> buffer(4096);
            while (remaining > 0) {
                auto batch = std::min(remaining,buffer.size());
                is.read(&buffer[0], batch);
                remaining -= batch;
                os.write(&buffer[0],batch);
            }
        }
        
        // copy the other file
        {
            std::ifstream  is(chunks.back()->filename());
            os << is.rdbuf();
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
        parse_size(options.buffer_size.getValue()),
        options.threads.getValue(),
        options.multiway_merge.getValue()
    };
    
}




    
