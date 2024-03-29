/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef ETAPDF
#define ETAPDF

#include "include/RooAbsPdf.h"
#include "include/RooRealProxy.h"
#include "include/RooCategoryProxy.h"
#include "include/RooAbsReal.h"
#include "include/RooAbsCategory.h"
 
class EtaPdf : public RooAbsPdf {
public:
  EtaPdf() {} ; 
  EtaPdf(const char *name, const char *title,
	      RooAbsReal& _theta);
  EtaPdf(const EtaPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new EtaPdf(*this,newname); }
  inline virtual ~EtaPdf() { }

protected:

  RooRealProxy theta ;
  
  Double_t evaluate() const ;

private:

  ClassDef(EtaPdf,1) // Your description goes here...
};
 
#endif
