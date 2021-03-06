// $Id$
//==============================================================================
//!
//! \file main_PoroElasticity.C
//!
//! \date April 16 2015
//!
//! \author Yared Bekele
//!
//! \brief Main program for the isogeometric solver for Poroelasticity
//!
//==============================================================================

#include "IFEM.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMPoroElasticity.h"
#include "SIMSolver.h"
#include "ASMmxBase.h"
#include "Utilities.h"
#include "Profiler.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "AppCommon.h"


  template<class Dim>
int runSimulator(char* infile, bool mixed)
{
  std::vector<unsigned char> fields;
  if (mixed)
    fields = {Dim::dimension, 1};
  else
    fields = {Dim::dimension+1};
  SIMPoroElasticity<Dim> model(fields);
  SIMSolver< SIMPoroElasticity<Dim> > solver(model);

  int res = ConfigureSIM(model,infile,true);
  if (res)
    return res;

  if (!solver.read(infile))
    return 3;

  // HDF5 output
  DataExporter* exporter = nullptr;
  if (model.opt.dumpHDF5(infile))
    exporter = SIM::handleDataOutput(model, solver, model.opt.hdf5,
                                     false, model.getDumpInterval(), 1);

  if (solver.solveProblem(infile, exporter))
    return 5;

  delete exporter;

  return 0;
}


int main(int argc, char ** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  std::vector<int> ignoredPatches;
  int i;
  char ndim = 3;
  char* infile = 0;
  bool mixed = false;

  IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      ndim = 2;
    else if (!strcmp(argv[i],"-mixed")) {
      ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
      mixed = true;
    }
    else if (!infile)
      infile = argv[i];
    else
      std::cerr << "*** Unknown option ignored: " << argv[i] << std::endl;

  if (!infile)
  {
    IFEM::cout << "Usage: " << argv[0]
               << " <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
               << " [-free] [-lag|-spec|-LR] [-1D|-2D] [-mixed] [-nGauss <n>]"
               << "\n       [-vtf <format> [-nviz <nviz>]"
               << " [-nu <nu> [-nv <nv>] [-nw <nw>]] [-hdf5]"<< std::endl;
    return 0;
  }

  IFEM::cout << "\n >>> IFEM Poroelasticity Solver <<<"
             << "\n =================================="
             << "\n Executing command:\n";
  for (i = 0; i < argc; i++) IFEM::cout << " " << argv[i];
  IFEM::cout << "\n\n Input file: " << infile;
  IFEM::getOptions().print(IFEM::cout);
  IFEM::cout << std::endl;

  if (ndim == 3)
    return runSimulator<SIM3D>(infile, mixed);
  else if (ndim == 2)
    return runSimulator<SIM2D>(infile, mixed);

  return 1;
}
