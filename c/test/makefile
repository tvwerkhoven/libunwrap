#
# This is a makefile to compile a stand-alone windows program.
#

ALL : unwraptest

INCLUDE=\
        -I"C:\Program Files\Microsoft Visual Studio 10.0\VC\include"  \
        -I"C:\Program Files\Microsoft SDKs\Windows\v7.0A\Include" \
        -I"..\unwrap"

#LIB=-L"C:\Program Files\Microsoft SDKs\Windows\v7.0A\Lib\" \
#       -L"C:\Program Files\Microsoft Visual Studio 10.0\VC\Lib\"


EXTRALIBS="C:\Program Files\Microsoft SDKs\Windows\v7.0A\Lib\Gdi32.lib"

SRC=..\unwrap\libunwrap.cc

CLEAN :
        -@erase unwrap-test.exe


unwraptest :
        cl.exe $(INCLUDE) $(SRC) test_unwrap.cc /EHsc /MDd /Zi /Feunwrap-test.exe

# minterf.lib 

#minterf :
#        cl.exe $(INCLUDE) $(DMLIBS) $(SHLIBS) $(IMGLIBS) $(MINTERSRC) $(EXTRALIBS) /LD /EHsc /MD

