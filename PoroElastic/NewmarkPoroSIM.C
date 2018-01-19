// $Id$
//==============================================================================
//!
//! \file NewmarkPoroSIM.C
//!
//! \date Jan 18 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Newmark solution driver for isogeometric dynamic poroelastic simulators.
//!
//==============================================================================

#include "NewmarkPoroSIM.h"
#include "SAM.h"
#include "SIMoutput.h"


bool NewmarkPoroSIM::solutionNorms (const TimeDomain&,
                                    double zero_tolerance, std::streamsize outPrec)
{
  if (msgLevel < 0 || solution.size() < 3) return true;

  // Cannot use the enums here because this method is inherited
  size_t a = solution.size()-1;
  size_t v = solution.size()-2;

  size_t d, nf = model.getNoFields(1);
  size_t iMax[nf], jMax[nf], kMax[nf], lMax[3];
  double dMax[nf], vMax[nf], aMax[nf], pMax[3];
  double disL2 = model.solutionNorms(solution[0],dMax,iMax,nf);
  double velL2 = model.solutionNorms(solution[v],vMax,jMax,nf);
  double accL2 = model.solutionNorms(solution[a],aMax,kMax,nf);
  double preL2, prrL2, praL2;
  if (model.mixedProblem()) {
    preL2 = model.solutionNorms(solution[0],pMax,  lMax,1,  'P');
    prrL2 = model.solutionNorms(solution[v],pMax+1,lMax+1,1,'P');
    praL2 = model.solutionNorms(solution[a],pMax+2,lMax+2,1,'P');
  } else {
    size_t len, len2;
    model.getSAM()->normL2(solution[0],'D',nf,len);
    double tmp = model.getSAM()->normL2(solution[0],'D',nf-1,len2);
    preL2 = sqrt((disL2*disL2*len-tmp*tmp*len2) / (len-len2));
    disL2 = tmp;
    tmp = model.getSAM()->normL2(solution[v],'D',nf-1,len2);
    prrL2 = sqrt((velL2*velL2*len-tmp*tmp*len2) / (len-len2));
    velL2 = tmp;
    tmp = model.getSAM()->normL2(solution[a],'D',nf-1,len2);
    praL2 = sqrt((accL2*accL2*len-tmp*tmp*len2) / (len-len2));
    accL2 = tmp;

    lMax[0] = iMax[nf-1];
    pMax[0] = dMax[nf-1];
    pMax[1] = vMax[nf-1];
    lMax[1] = jMax[nf-1];
    pMax[2] = aMax[nf-1];
    lMax[2] = kMax[nf-1];
    --nf;
  }

  utl::LogStream& cout = model.getProcessAdm().cout;
  std::streamsize stdPrec = outPrec > 0 ? cout.precision(outPrec) : 0;
  double old_tol = utl::zero_print_tol;
  utl::zero_print_tol = zero_tolerance;
  char D;

  cout <<"  Displacement L2-norm            : " << utl::trunc(disL2);
  for (d = 0, D = 'X'; d < nf; d++, D=='Z' ? D='x' : D++)
    if (utl::trunc(dMax[d]) != 0.0)
      cout <<"\n               Max "<< D
           <<"-displacement : "<< dMax[d] <<" node "<< iMax[d];
  cout << "\n  Pressure L2-norm                : " << utl::trunc(preL2);
  if (utl::trunc(pMax[0]) != 0.0)
      cout <<"\n               Max pressure       : "
           <<pMax[0] <<" node "<< lMax[0];

  cout <<"\n  Velocity L2-norm                : "<< utl::trunc(velL2);
  for (d = 0, D = 'X'; d < nf; d++, D=='Z' ? D='x' : D++)
    if (utl::trunc(vMax[d]) != 0.0)
      cout <<"\n               Max "<< D
           <<"-velocity     : "<< vMax[d] <<" node "<< jMax[d];
  cout << "\n  Pressure rate L2-norm           : " << utl::trunc(prrL2);
  if (utl::trunc(pMax[1]) != 0.0)
      cout <<"\n               Max pressure rate  : "
           <<pMax[1] <<" node "<< lMax[1];

  cout <<"\n  Acceleration L2-norm            : "<< utl::trunc(accL2);
  for (d = 0, D = 'X'; d < nf; d++, D=='Z' ? D='x' : D++)
    if (utl::trunc(aMax[d]) != 0.0)
      cout <<"\n               Max "<< D
           <<"-acceleration : "<< aMax[d] <<" node "<< kMax[d];
  cout << "\n  Pressure acceleration L2-norm   : " << utl::trunc(praL2);
  if (utl::trunc(pMax[2]) != 0.0)
      cout <<"\n        Max pressure acceleration : "
           <<pMax[2] <<" node "<< lMax[2];

  cout << std::endl;
  utl::zero_print_tol = old_tol;
  if (stdPrec > 0) cout.precision(stdPrec);

  return true;
}
