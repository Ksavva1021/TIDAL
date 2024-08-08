#ifndef MULTIDRAW_H
#define MULTIDRAW_H

#include "Rtypes.h"

class TTree;
class TTreeFormula;
class TObjArray;

void MultiDraw(TTree *inTree,
               TObjArray *Formulae, TObjArray *Weights, TObjArray *Hists,
               UInt_t ListLen);

#endif // MULTIDRAW_H
