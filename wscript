#! /usr/bin/env python
# encoding: utf-8

def options(opt):
    opt.load('compiler_cxx')
    opt.load('compiler_c')


#    ctx.add_option('--tcmalloc', action='store_true', default=False,
#                help='execute the program after it is built')

def configure(cnf):
    cnf.load('compiler_cxx')
    cnf.load('compiler_c')
    cnf.check(features='cxx cxxprogram',msg='cxx and cxxprogram') 

def build(bld):
    bld.recurse('src')        # programs





#     bld.recurse('tessellate') # libtessellate




    # conf.env.CFLAGS = ['-g']
    # conf.setenv('release')
    # conf.load('compiler_c')
    # conf.env.CFLAGS = ['-O2']



#, lib=['m'], cflags=['-Wall'], defines=['var=foo'], uselib_store='M')

# def build(bld):
#     # bld(features='c cshlib', source='b.c', target='mylib')
#     # bld(features='cxx cxxprogram', source='rastahit.cc', target='app') #, use=['M','mylib'], lib=['dl'])
#     bld.program(source="rastahit.cc", target='app')

