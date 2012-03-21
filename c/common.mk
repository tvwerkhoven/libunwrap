## @file c/common.mk
## @author Tim van Werkhoven (werkhoven@strw.leidenuniv.nl)
## Common makefile directives for all subdirs.

### Inclusion directories

LIB_DIR = $(top_builddir)/unwrap

### Common flags

AM_CPPFLAGS = -I$(LIB_DIR)
#				$(COMMON_CFLAGS)

AM_CPPFLAGS += -D__STDC_FORMAT_MACROS \
		-D__STDC_LIMIT_MACROS
    
#LDADD = $(COMMON_LIBS)

# More error reporting during compilation
AM_CFLAGS = -Wall -Wextra -Wfatal-errors

### Debug options
if HAVE_DEBUG
AM_CFLAGS += -ggdb -g3 -O0 -fno-inline
else
AM_CFLAGS += -O3 -ftree-vectorize
endif
