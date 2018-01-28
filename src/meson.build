# File organization
libcrcfastSrc   = files('crcfast.c')
crcfastmainSrc  = files('crcfastmain.c')

#
# Get C compiler
#
cc = meson.get_compiler('c')


#
# Compiler args
#
libcrcfastCArgs = []
if cc.has_argument('-mssse3')
    libcrcfastCArgs += ['-mssse3']
endif


#
# Dynamic Libraries need careful consideration of their compiler arguments!
# Static libraries, not so much.
#
libcrcfastShaCArgs  = libcrcfastCArgs
libcrcfastShaCArgs += ['-DCRCFAST_SHAREDOBJECT']
if host_machine.system() == 'windows'
    libcrcfastShaCArgs += ['-DSICH_BUILDING_DLL']
else
    libcrcfastShaCArgs += ['-fvisibility=hidden']
endif

libcrcfastStaCArgs  = libcrcfastCArgs
libcrcfastStaCArgs += []



# Build libraries
libcrcfastShaInstall = get_option('default_library') == 'shared'
libcrcfastSha = shared_library('crcfast',           libcrcfastSrc,
                               include_directories: libcrcfastInc,
                               c_args:              libcrcfastShaCArgs,
                               link_args:           [],
                               install:             libcrcfastShaInstall,)

libcrcfastStaInstall = get_option('default_library') == 'static'
libcrcfastSta = static_library('crcfast',           libcrcfastSrc,
                               include_directories: libcrcfastInc,
                               c_args:              libcrcfastStaCArgs,
                               link_args:           [],
                               install:             libcrcfastStaInstall,)


# Build executable
crcfastmainExe = executable   ('crcfastmain',       crcfastmainSrc,
                               include_directories: libcrcfastInc,
                               c_args:              [],
                               link_whole:          [libcrcfastSta],
                               install:             false,)

