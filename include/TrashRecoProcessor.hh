#ifndef TrashRecoProcessor_h
#define TrashRecoProcessor_h 1
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <EVENT/LCCollection.h>
#include <EVENT/Vertex.h>
#include <EVENT/ReconstructedParticle.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "marlin/Processor.h"
#include "lcio.h"

#include "MathOperator.hh"
#include "VertexRecoOperator.hh"
#include "JetOperator.hh"
#include "VertexTag.hh"
#include <string>
#include <vector>
#include <TFile.h>
//#include <TLorentzVector.h>
#include <TTree.h>

using namespace lcio ;
using namespace marlin ;


namespace TTbarAnalysis 
{
	class TrashRecoProcessor : public Processor {
	  
	 public:
	  
	  virtual Processor*  newProcessor() { return new TrashRecoProcessor ; }
	  
	  
	  TrashRecoProcessor() ;
	  
	  /** Called at the begin of the job before anything is read.
	   * Use to initialize the processor, e.g. book histograms.
	   */
	  virtual void init() ;
	  
	  /** Called for every run.
	   */
	  virtual void processRunHeader( LCRunHeader* run ) ;
	  
	  /** Called for every event - the working horse.
	   */
	  virtual void processEvent( LCEvent * evt ) ; 
	  
	  
	  virtual void check( LCEvent * evt ) ; 
	  
	  
	  /** Called after data processing for clean up.
	   */
	  virtual void end() ;
	  void Write(VertexTag * tag, int number);
	  void Write(Vertex * vertex, int number);
	  void Write(std::vector< Jet * > * jets);
	  Vertex * FindPrimaryVertex(const LCCollection * collection); 
	  void PrintParticle(ReconstructedParticle * particle);
	  void WriteTagged(std::vector< VertexTag * > * tagged);
	  float getMissingPt(std::vector< ReconstructedParticle * > & bdaugthers, std::vector< VertexTag * > * tags, int pdg);
	  VertexTag * getBvertex(std::vector< VertexTag * > * tags, int pdg);
	  void ClearVariables();
	  void GetBAngles(std::vector< VertexTag * > * tags);
	 protected:
	
	  /** Input collection name.
	   */
	  std::string _colSecName ;
	  std::string _colPriName ;
	  std::string _colMCName;
	  std::string _colJetsName;
	  std::string _colJetsRelName;
	  std::string _colquarkName;
	  
	  int _nRun ;
	  int _nEvt ;
	
	  TFile * _hfile;
	  TTree * _hTree;
	  TTree * _hTaggedTree;
	  TTree * _hUntaggedTree;
	  TTree * _hJetTree;
	  std::string _hfilename ;
	  float _angleAcceptance;	  
	  int _handleJets; 

	  int _numberOfTagged;
	  int _numberOfTotal;
	  int _numberOfTernary;
	  int _numberOfSecondary;
	  int _numberOfUnknown;
	  int _bnumber1;
	  int _bbarnumber1;
	  float _bmomentum;
	  float _bbarmomentum;
	  int _bexists;
	  int _bbarexists;
	  int _bcharge;
	  int _bbarcharge;
	  float _bptmiss;
	  float _bbarptmiss;
	  float _bbarteta;
	  float _bteta;
	  float _bIPdistance;
	  float _bbarIPdistance;
	  float _bbarchimean;
	  float _bchimean;
	  float _bprobmean;
	  float _bbarprobmean;
	  static const int MAXV2 = 36;
	  float _distances[MAXV2];
	  int _numberOfDistances;
	  static const int MAXV = 15;
	  float _distanceFromIP[MAXV];
	  float _coordinates[MAXV][3];
	  int _PDG[MAXV];
	  float _probability[MAXV];
	  int _charge[MAXV];
	  float _mindistance[MAXV];
	  float _chi2[MAXV];
	  int _generation[MAXV];
	  int _numberOfParticles[MAXV];
	  float _energyOfParticles[MAXV][MAXV];
	  float _momentumOfParticles[MAXV][MAXV];
	  float _massOfParticles[MAXV][MAXV];
	  //TLorentzVector * _particles[MAXV][MAXV]
	  Vertex * _primary;
	  
	  int _numberOfJets;
	  float _btags[MAXV];
	  float _ctags[MAXV];
	  int _mcpdg[MAXV];
	  int _nvertices[MAXV];

	
	} ;
		
} /* TTbarAnalysis */

#endif



