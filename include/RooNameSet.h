/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooNameSet.h,v 1.16 2007/05/11 09:11:30 verkerke Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_NAME_SET
#define ROO_NAME_SET

#include "TString.h"
#include "RooAbsArg.h"
#include "RooPrintable.h"
class RooArgSet ;

class RooNameSet : public TObject, public RooPrintable {
public:

  // Constructors, assignment etc.
  RooNameSet();
  RooNameSet(const RooArgSet& argSet);
  RooNameSet(const RooNameSet& other) ;
  virtual TObject* Clone(const char*) const { return new RooNameSet(*this) ; }
  virtual ~RooNameSet() ;

  void refill(const RooArgSet& argSet) ;
  RooArgSet* select(const RooArgSet& list) const ;
  Bool_t operator==(const RooNameSet& other) ;  
  RooNameSet& operator=(const RooNameSet&) ;

  virtual void printName(ostream& os) const ;
  virtual void printTitle(ostream& os) const ;
  virtual void printClassName(ostream& os) const ;
  virtual void printValue(ostream& os) const ;

  inline virtual void Print(Option_t *options= 0) const {
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

  void setNameList(const char* givenList) ;

  Int_t _len ;
  char* _nameList ; //[_len]

protected:

  void extendBuffer(Int_t inc) ;

  ClassDef(RooNameSet,1) // A sterile version of RooArgSet, containing only the names of the contained RooAbsArgs
};

#endif
