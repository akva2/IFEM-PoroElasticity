// $Id$
//==============================================================================
//!
//! \file PoroElasticity.C
//!
//! \date April 16 2015
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for time-dependent poroelasticity problems.
//!
//==============================================================================

#include "PoroElasticity.h"
#include "PoroMaterial.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml.h"


//! \brief Enum for element level solution vectors
enum SolutionVectors
{
  U = 0, // displacement
  P = 1, // pore pressure
  NSOL = 2
};


//! \brief Enum for element level right-hand-side vectors
enum ResidualVectors
{
  Fu = 0,
  Fp = 1,
  Fres = 2,
  Fprev = 3,
  NVEC = 4
};


//! \brief Enum for element level left-hand-side matrices
enum TangentMatrices
{
  uu = 0,                       // Stiffness matrix
  up = 1,                       // Coupling matrix
  pp_S = 2,                     // Compressibility matrix
  pp_P = 3,                     // Permeability matrix
  pp = 4,                       // S + dt P
  Ktan = 5,
  Kprev = 6,
  NMAT = 7
};


PoroElasticity::MixedElmMats::MixedElmMats ()
{
  this->resize(NMAT,NVEC); // Number of element matrices and vectors
}


void PoroElasticity::MixedElmMats::makeNewtonMatrix_P(Matrix& N, size_t pp_idx) const
{
  size_t n = A[uu].rows();
  N.fillBlock(A[up], 1+n, 1, true);
  N.fillBlock(A[pp_idx], 1+n, 1+n);
}


void PoroElasticity::MixedElmMats::makeNewtonMatrix_U(Matrix& N) const
{
  size_t n = A[uu].rows();
  N.fillBlock(A[uu], 1, 1);
  N.addBlock(A[up], -1.0, 1, 1+n);
}


const Matrix& PoroElasticity::MixedElmMats::getNewtonMatrix () const
{
  Matrix& N = const_cast<Matrix&>(A[Ktan]);

  makeNewtonMatrix_U(N);
  makeNewtonMatrix_P(N, pp);

  return A[Ktan];
}


const Vector& PoroElasticity::MixedElmMats::getRHSVector () const
{
  Vector& F = const_cast<Vector&>(b[Fres]);
  size_t n = b[Fu].size();

  std::copy(b[Fu].begin(), b[Fu].end(), F.begin());
  std::copy(b[Fp].begin(), b[Fp].end(), F.begin()+n);

  F += b[Fprev];

  return b[Fres];
}


PoroElasticity::NonMixedElmMats::NonMixedElmMats ()
{
  this->resize(NMAT,NVEC); // Number of element matrices and vectors
}


void PoroElasticity::NonMixedElmMats::makeNewtonMatrix_P(Matrix& N, size_t pp_idx) const
{
  size_t n = A[pp].rows();
  size_t nsd = A[uu].rows()/A[pp].rows();
  size_t nf = nsd + 1;

  for (size_t i = 1; i <= n; ++i) {
    for (size_t j = 1; j <= n; ++j) {
      // PP-coupling
      N(nf*i, nf*j) = A[pp_idx](i,j);
      for (size_t l = 1; l <= nsd; ++l) {
        // PU-coupling
        N(j*nf, nf*(i-1)+l) = A[up](nsd*(i-1)+l, j);
      }
    }
  }
}


void PoroElasticity::NonMixedElmMats::makeNewtonMatrix_U(Matrix& N) const
{
  size_t n = A[pp].rows();
  size_t nsd = A[uu].rows()/A[pp].rows();
  size_t nf = nsd + 1;

  for (size_t i = 1; i <= n; ++i) {
    for (size_t j = 1; j <= n; ++j) {
      for (size_t l = 1; l <= nsd; ++l) {
        // UP-coupling
        N(nf*(i-1)+l, j*nf) = -A[up](nsd*(i-1)+l, j);
        for (size_t k = 1; k <= nsd; ++k) {
          // UU-coupling
          N(nf*(i-1)+l, (j-1)*nf+k) = A[uu](nsd*(i-1)+l, (j-1)*nsd+k);
        }
      }
    }
  }
}


const Matrix& PoroElasticity::NonMixedElmMats::getNewtonMatrix () const
{
  Matrix& N = const_cast<Matrix&>(A[Ktan]);

  makeNewtonMatrix_U(N);
  makeNewtonMatrix_P(N, pp);

  return A[Ktan];
}


const Vector& PoroElasticity::NonMixedElmMats::getRHSVector () const
{
  Vector& F = const_cast<Vector&>(b[Fres]);

  size_t nsd = b[Fu].size()/b[Fp].size();
  utl::interleave(b[Fu], b[Fp], F, nsd, 1);

  F += b[Fprev];

  return b[Fres];
}


PoroElasticity::PoroElasticity (unsigned short int n) : Elasticity(n)
{
  sc = 0.0;
  gacc = 9.81; //kmo: Danger! hard-coded physical property. Why not derive this one from gravity.length() instead ???
}


bool PoroElasticity::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"scaling"))
    return this->Elasticity::parse(elem);
  else if (utl::getAttribute(elem,"value",sc))
    IFEM::cout <<"\tScaling: sc = "<< sc << std::endl;

  return true;
}


Material* PoroElasticity::parseMatProp (const TiXmlElement* elem, bool)
{
  IFEM::cout <<" Poroelastic material, see \"Problem definition:\" below."
             << std::endl;
  material = new PoroMaterial();
  material->parse(elem);
  return material;
}


void PoroElasticity::printLog () const
{
  IFEM::cout <<"PoroElasticity: scaling = "<< sc << std::endl;
  this->Elasticity::printLog();
}


void PoroElasticity::initLocalIntegral(ElmMats *result, size_t ndof_displ,
                                       size_t ndof_press, bool neumann) const
{
  size_t ndof_tot = ndof_displ + ndof_press;

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->b[Fres].resize(ndof_tot);
  result->b[Fprev].resize(ndof_tot);
  result->b[Fu].resize(ndof_displ);
  result->b[Fp].resize( ndof_press);

  if (!neumann)
  {
    result->A[uu].resize(ndof_displ, ndof_displ);
    result->A[up].resize(ndof_displ, ndof_press);
    result->A[pp_S].resize(ndof_press, ndof_press);
    result->A[pp_P].resize(ndof_press, ndof_press);
    result->A[pp].resize(ndof_press, ndof_press);
    result->A[Ktan].resize(ndof_tot, ndof_tot);
    result->A[Kprev].resize(ndof_tot, ndof_tot);
  }
}


LocalIntegral* PoroElasticity::getLocalIntegral (const std::vector<size_t>& nen,
                                                 size_t, bool neumann) const
{
  ElmMats* result = new MixedElmMats();
  initLocalIntegral(result, nsd * nen[0], nen[1], neumann);
  return result;
}


LocalIntegral* PoroElasticity::getLocalIntegral (size_t nen,
                                                 size_t, bool neumann) const
{
  ElmMats* result = new NonMixedElmMats();
  initLocalIntegral(result, nsd * nen, nen, neumann);
  return result;
}


bool PoroElasticity::initElement (const std::vector<int>& MNPC,
                                  const std::vector<size_t>& elem_sizes,
                                  const std::vector<size_t>& basis_sizes,
                                  LocalIntegral& elmInt)
{
  if (primsol.empty() || primsol.front().empty())
    return true;

  // Split the nodal correspondance array into one for each solution field
  std::vector<int>::const_iterator pStart = MNPC.begin() + elem_sizes.front();
  std::vector<int> MNPCu(MNPC.begin(),pStart);
  std::vector<int> MNPCp(pStart,MNPC.end());

  // Extract the element level solution vectors
  if (elmInt.vec.size() < NSOL) elmInt.vec.resize(NSOL);
  int ierr = utl::gather(MNPCu, nsd, primsol.front(), elmInt.vec[U])
           + utl::gather(MNPCp, 0,1, primsol.front(), elmInt.vec[P],
                         nsd*basis_sizes.front(), basis_sizes.front());
  if (ierr == 0) return true;

  std::cerr <<" *** PoroElasticity::initElement: Detected "<< ierr/2
            <<" node numbers out of range."<< std::endl;

  return false;
}


bool PoroElasticity::initElement (const std::vector<int>& MNPC,
                                  LocalIntegral& elmInt)
{
  if (primsol.empty() || primsol.front().empty())
    return true;

  // Extract the element level solution vectors
  if (elmInt.vec.size() < NSOL) elmInt.vec.resize(NSOL);
  Matrix temp(nsd+1,MNPC.size());
  int ierr = utl::gather(MNPC, nsd+1, primsol.front(), temp);
  if (ierr == 0)
  {
    elmInt.vec[P] = temp.getRow(nsd+1);
    elmInt.vec[U] = temp.expandRows(-1);
    return true;
  }

  std::cerr <<" *** PoroElasticity::initElement: Detected "<< ierr
            <<" node numbers out of range."<< std::endl;

  return false;
}


bool PoroElasticity::evalIntMx (LocalIntegral& elmInt,
                                const MxFiniteElement& fe,
                                const TimeDomain& time, const Vec3& X) const
{
  return this->evalInt(elmInt, fe, time, X);
}


bool PoroElasticity::evalElasticityMatrices (ElmMats& elMat, const Matrix& B,
                                             const FiniteElement& fe,
                                             const Vec3& X) const
{
  Matrix C;
  SymmTensor eps(nsd), sigma(nsd); double U = 0.0;
  if (!material->evaluate(C,sigma,U,fe,X,eps,eps,0))
    return false;

  Matrix CB;
  CB.multiply(C,B).multiply(fe.detJxW); // CB = dSdE*B*|J|*w
  elMat.A[uu].multiply(B,CB,true,false,true); // EK += B^T * CB

  return true;
}


bool PoroElasticity::evalCouplingMatrix (Matrix& mx, const Matrix& B,
                                         const Vector& N, double scl) const
{
  Matrix K(N.size(), nsd * (nsd + 1) / 2);
  for (size_t i = 1; i <= N.size(); i++)
    for (size_t j = 1; j <= nsd; j++)
      K(i,j) = scl * N(i);

  mx.multiply(B, K, true, true, true);

  return true;
}


bool PoroElasticity::evalCompressibilityMatrix (Matrix& mx, const Vector& N,
                                                double scl) const
{
  Matrix temp(mx.rows(), mx.cols());
  if (!temp.outer_product(N,N))
    return false;

  mx.add(temp,scl);

  return true;
}


bool PoroElasticity::evalPermeabilityMatrix (Matrix& mx, const Matrix& dNdX,
                                             const Vec3& permeability,
                                             double scl) const
{
  for (size_t i = 1; i <= dNdX.rows(); i++)
    for (size_t j = 1; j <= dNdX.rows(); j++)
      for (size_t k = 1; k <= nsd; k++)
        mx(i,j) += scl * permeability[k-1]  * dNdX(i,k) * dNdX(j,k);

  return true;
}


bool PoroElasticity::evalInt (LocalIntegral& elmInt,
                              const FiniteElement& fe,
                              const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  const PoroMaterial* pmat = dynamic_cast<const PoroMaterial*>(material);
  if (!pmat) {
    std::cerr << __FUNCTION__ << ": No material data." << std::endl;
    return false;
  }

  Matrix Bmat;
  if (!this->formBmatrix(Bmat,fe.dNdX))
    return false;

  Vec3 permeability = pmat->getPermeability(X);

  double rhog = pmat->getFluidDensity(X) * gacc;
  double scl = sc;
  if (scl == 0.0)
    scl = sqrt(pmat->getStiffness(X)*rhog/(permeability.x*time.dt));

  // Biot's coefficient
  double Ko = pmat->getBulkMedium(X);
  double Ks = pmat->getBulkSolid(X);
  double Kw = pmat->getBulkWater(X);
  double poro = pmat->getPorosity(X);

  double alpha = 1.0 - (Ko/Ks);
  // Inverse of the compressibility modulus
  double Minv = ((alpha - poro)/Ks) + (poro/Kw);

  if (!this->evalElasticityMatrices(elMat, Bmat, fe, X))
    return false;

  if (!this->evalCouplingMatrix(elMat.A[up], Bmat, fe.basis(2),
                                scl*alpha*fe.detJxW))
    return false;

  if (!this->evalCompressibilityMatrix(elMat.A[pp_S], fe.basis(2), 
                                       scl*scl*Minv*fe.detJxW))
    return false;

  return this->evalPermeabilityMatrix(elMat.A[pp_P], fe.grad(2),
                                      permeability, scl*scl/rhog*fe.detJxW);
}


bool PoroElasticity::evalBouMx (LocalIntegral& elmInt,
                                const MxFiniteElement& fe,
                                const TimeDomain&, const Vec3& X,
                                const Vec3& normal) const
{
  // Using the inherited standard method from Elasticity
  return this->evalBou(elmInt, fe, X, normal);
}


bool PoroElasticity::finalizeElement (LocalIntegral& elmInt, const TimeDomain& time, size_t)
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  elMat.A[pp].add(elMat.A[pp_S]);
  elMat.A[pp].add(elMat.A[pp_P], time.dt);

  Vector prevSol;

  MixedElmMats* mMat = dynamic_cast<MixedElmMats*>(&elmInt);
  if (mMat) {
    mMat->makeNewtonMatrix_P(elMat.A[Kprev], pp_S);
    prevSol = elmInt.vec[U];
    prevSol.insert(prevSol.end(), elmInt.vec[P].begin(), elmInt.vec[P].end());
  } else {
    NonMixedElmMats* mMat = dynamic_cast<NonMixedElmMats*>(&elmInt);
    mMat->makeNewtonMatrix_P(elMat.A[Kprev], pp_S);
    utl::interleave(elmInt.vec[U], elmInt.vec[P], prevSol, nsd, 1);
  }

  elMat.b[Fprev] = elMat.A[Kprev] * prevSol;

  return true;
}


size_t PoroElasticity::getNoFields (int fld) const
{
  if (fld < 2)
    return nsd+1;
  return nsd * (nsd + 1);
}


std::string PoroElasticity::getField1Name (size_t i, const char* prefix) const
{
  if (i == 11)
    return "Displacements";
  if (i == 12)
    return "Pressure";

  if (i >= nsd)
    i = 3;

  static const char* s[5] = {"u_x", "u_y", "u_z", "p^w"};

  if (!prefix)
    return s[i];
  return prefix + std::string(" ") + s[i];
}


std::string PoroElasticity::getField2Name (size_t i, const char* prefix) const
{
  static const char* s[][6] = {{"x", "y", "xy"},
                               {"x", "y", "z", "yz", "xz", "xy"}};
  size_t ncomps = nsd * (nsd + 1) / 2;

  std::string name = (i < ncomps ? "eps" : "sig") + std::string("_") + s[nsd-2][i % ncomps];
  if (!prefix)
    return name;

  return prefix + std::string(" ") + name;
}


bool PoroElasticity::evalSol (Vector& s, const MxFiniteElement& fe,
                              const Vec3& X, const std::vector<int>& MNPC,
                              const std::vector<size_t>& elem_sizes) const
{
  Vector eV;
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes.front();
  utl::gather(std::vector<int>(MNPC.begin(),fstart),nsd,primsol.front(),eV);

  return this->evalSolCommon(s,fe,X,Vector(eV.ptr(),nsd*fe.N.size()));
}


bool PoroElasticity::evalSol (Vector& s, const FiniteElement& fe,
                              const Vec3& X, const std::vector<int>& MNPC) const
{
  Vector eV;
  utl::gather(MNPC,nsd+1,primsol.front(),eV);

  Vector disp(nsd * fe.N.size());
  for (size_t i = 0; i < nsd; i++)
    for (size_t bfun = 0; bfun < fe.N.size(); bfun++)
      disp[nsd*bfun+i] = eV[(nsd+1)*bfun+i];

  return this->evalSolCommon(s,fe,X,disp);
}


bool PoroElasticity::evalSolCommon (Vector& s,
                                    const FiniteElement& fe, const Vec3& X,
                                    const Vector& disp) const
{
  if (!material)
  {
    std::cerr << __FUNCTION__ <<": No material data."<< std::endl;
    return false;
  }

  Matrix Bmat;
  if (!this->formBmatrix(Bmat,fe.dNdX))
    return false;

  SymmTensor eps(nsd), sigma(nsd);
  if (!Bmat.multiply(disp,eps))
    return false;

  Matrix Cmat; double U = 0.0;
  if (!material->evaluate(Cmat,sigma,U,fe,X,eps,eps))
    return false;

  s = eps;
  const RealArray& sig = sigma;
  s.insert(s.end(),sig.begin(),sig.end());

  return true;
}
