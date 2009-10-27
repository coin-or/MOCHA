# Copyright 2009 by Jesus De Loera, David Haws, Jon Lee, Allison
# O'Hair, University of California at Davis, IBM, Inc.
# All Rights Reserved.
# This file is distributed under the Eclips Public License v 1.0

## $Id$

# Authors:  David Haws      University of California at Davis   2009-10-27

AUTOMAKE_OPTIONS = foreign

bin_PROGRAMS = genrandmo matroidtest graphtest testmatrix mvbalclust nagibatest genmatrix alltoproj point2po designtomatrix estimatebases tabusearch localsearch

genrandmo_SOURCES = genrandmo.cpp mathprog.cpp matrix.cpp graph.cpp matroid.cpp
genrandmo_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#genrandmo_LDFLAGS = -lblas

matroidtest_SOURCES = matroidtest.cpp matrix.cpp graph.cpp matroid.cpp mathprog.cpp
matroidtest_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#matroidtest_LDFLAGS = -lblas

graphtest_SOURCES = graphtest.cpp matrix.cpp graph.cpp
graphtest_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#graphtest_LDFLAGS = -lblas

testmatrix_SOURCES = testmatrix.cpp matrix.cpp
testmatrix_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#testmatrix_LDFLAGS = -lblas

mvbalclust_SOURCES = mvbalclust.cpp mathprog.cpp matrix.cpp graph.cpp matroid.cpp 
mvbalclust_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#mvbalclust_LDFLAGS = -lblas

nagibatest_SOURCES = nagibatest.cpp graph.cpp matrix.cpp
nagibatest_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#nagibatest_LDFLAGS = -lblas

genmatrix_SOURCES = genmatrix.cpp matrix.cpp
genmatrix_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#genmatrix_LDFLAGS = -lblas

alltoproj_SOURCES = alltoproj.cpp matrix.cpp
alltoproj_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#alltoproj_LDFLAGS = -lblas

point2po_SOURCES = point2po.cpp mathprog.cpp matrix.cpp graph.cpp matroid.cpp 
point2po_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#point2po_LDFLAGS = -lblas

designtomatrix_SOURCES = designtomatrix.cpp matrix.cpp
designtomatrix_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#designtomatrix_LDFLAGS = -lblas

estimatebases_SOURCES = estimatebases.cpp matroid.cpp matrix.cpp graph.cpp
estimatebases_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#estimatebases_LDFLAGS = -lblas

tabusearch_SOURCES = tabusearch.cpp matroid.cpp matrix.cpp graph.cpp mathprog.cpp
tabusearch_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#tabusearch_LDFLAGS = -lblas

localsearch_SOURCES = localsearch.cpp matroid.cpp matrix.cpp graph.cpp mathprog.cpp
localsearch_LDADD = $(THIRDBUILDLIBS) $(ADDLIBS)
#localsearch_LDFLAGS = -lblas

EXTRA_DIST = matrix.h matroid.h graph.h mathprog.h
