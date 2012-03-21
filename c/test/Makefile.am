## @file test/Makefile.am
## @author Tim van Werkhoven (werkhoven@strw.leidenuniv.nl)

if WITH_TESTING

include $(top_srcdir)/common.mk

AM_CFLAGS += -I${top_srcdir}/unwrap/ -L${top_srcdir}/unwrap/

noinst_PROGRAMS = unwrap-test

unwrap_test_SOURCES = test_unwrap.c
unwrap_test_LDADD = ${top_srcdir}/unwrap/libunwrap.a \
    $(LDADD)

endif WITH_TESTING

