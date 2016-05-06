#ifndef TrashRecoProcessor_h
#define TrashRecoProcessor_h 1
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/Vertex.h>
#include <IMPL/VertexImpl.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "marlin/Processor.h"
#include "lcio.h"

#include "AlgebraImplementation.hh"
#include "MathOperator.hh"
#include "TrackOperator.hh"
#include "VertexRecoOperator.hh"
#include "JetOperator.hh"
#include "ParticleOperator.hh"
#include "JetVertexOperator.hh"
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
	  void Write(VertexTag * tag);
	  void Write(Vertex * vertex);
	  void WriteMissed(Vertex * vertex, LCCollection * collection);
	  void Write(std::vector< Jet * > * jets, LCCollection * rel, LCCollection * mc);
	  Vertex * FindPrimaryVertex(const LCCollection * collection); 
	  void PrintParticle(ReconstructedParticle * particle);
	  LCCollection * WriteTagged(std::vector< Jet * > * jets);
	  void WriteVertex(const std::vector< VertexTag * > & tags, LCCollectionVec * reco);
	  void WriteTagged(std::vector< VertexTag * > * tagged);
	  void WriteTaggedCollection(LCEvent * evt, std::vector< VertexTag * > * tags);
	  //float getMissingPt(std::vector< ReconstructedParticle * > & bdaugthers, std::vector< VertexTag * > * tags, int pdg);
	  float getMissingPt(const std::vector< VertexTag * > & tags, int pdg);
	  VertexTag * getBvertex(const std::vector< VertexTag * > & tags, int pdg);
	  int getTracksWithOffsets( Jet * jet);
	  void ClearVariables();
	  void AnalyseSecondaries(const LCCollection * prongs,const LCCollection * rel, const LCCollection * reco);
	  void Write(LCEvent * evt, std::vector< Particle > * missed, std::vector< ReconstructedParticle * > * m = NULL);
	  void GetBAngles(const std::vector< VertexTag * > & tags);
	  float GetDeltaP(ReconstructedParticle * particle);
	 protected:
	
	  /** Input collection name.
	   */
	  std::string _colSecName ;
	  std::string _colPriName ;
	  std::string _colMCName;
	  std::string _colPFOName;
	  std::string _colProngsName;
	  std::string _colRelName;
	  std::string _colTrackRelName;
	  std::string _colJetsName;
	  std::string _colJetsRelName;
	  std::string _colquarkName;
	  std::string _colMCMissName;
	  std::string _colMissName;
	  std::string _colMissVtxName;
	  std::string _colTagName;
	  std::string _colRecoProngsName;
	  std::string _colRecoProngsTracksName;
	  int _nRun ;
	  int _nEvt ;
	
	  TFile * _hfile;
	  TTree * _hTree;
	  TTree * _hTaggedTree;
	  TTree * _hSecTree;
	  TTree * _hUntaggedTree;
	  TTree * _hJetTree;
	  TTree * _hMissedTree;
	  TTree * _hMissedVertexTree;
	  std::string _hfilename ;
	  float _angleAcceptance;	  
	  int _handleJets; 
	  TrackOperator myTrackOperator;
	  int _numberOfTagged;
	  int _numberOfTotal;
	  int _numberOfTernary;
	  int _numberOfSecondary;
	  int _numberOfUnknown;
	  int _bnumber1;
	  int _bbarnumber1;
	  int _bgennumber;
	  int _bbargennumber;
	  int _bgencharge;
	  int _bbargencharge;
	  int _bnoffsettracks;
	  int _bbarnoffsettracks;
	  float _btag;
	  float _bbartag;
	  float _bmomentum;
	  float _bbarmomentum;
	  int _bnvtx;
	  int _bbarnvtx;
	  int _bcharge;
	  int _bbarcharge;
	  float _bptmiss;
	  float _bbarptmiss;
	  float _bbarcostheta;
	  float _bbarmass;
	  float _bmass;
	  float _bcostheta;
	  float _bIPdistance;
	  float _bbarIPdistance;
	  float _bbarchimean;
	  float _bchimean;
	  float _bprobmean;
	  float _bbarprobmean;
	  static const int MAXV2 = 36;
	  static const int MAXV3 = 100;
	  float _distances[MAXV2];
	  int _numberOfDistances;
	  static const int MAXV = 15;
	  float _distanceFromIP[MAXV];
	  float _precision[MAXV];
	  float _precisionT[MAXV];
	  float _distanceIP[MAXV];
	  float _coordinates[MAXV][3];
	  int _PDG[MAXV];
	  float _probability[MAXV];
	  float _angle[MAXV];
	  float _costhetaVtx[MAXV];
	  int _charge[MAXV];
	  int _truthNumber[MAXV];
	  float _mindistance[MAXV];
	  float _chi2[MAXV];
	  int _status[MAXV];
	  int _generation[MAXV];
	  int _numberOfParticles[MAXV];
	  float _energyOfParticles[MAXV][MAXV];
	  float _momentumOfParticles[MAXV][MAXV];
	  float _thetaOfParticles[MAXV][MAXV];
	  float _ptOfParticles[MAXV][MAXV];
	  float _costhetaOfParticles[MAXV][MAXV];
	  float _phiOfParticles[MAXV][MAXV];
	  int _typeOfParticles[MAXV][MAXV];
	  float _angleOfParticles[MAXV][MAXV];
	  float _offsetOfParticles[MAXV][MAXV];
	  float _deviationOfParticles[MAXV][MAXV];
	  float _chi2OfParticles[MAXV][MAXV];
	  float _errorOfParticles[MAXV][MAXV];
	  float _deltapOfParticles[MAXV][MAXV];
	  float _deltaStartOfParticles[MAXV][MAXV];
	  float _phi0OfParticles[MAXV][MAXV];
	  float _d0OfParticles[MAXV][MAXV];
	  float _z0OfParticles[MAXV][MAXV];
	  float _d0PrimeOfParticles[MAXV][MAXV];
	  float _rHitOfParticles[MAXV][MAXV];
	  float _z0PrimeOfParticles[MAXV][MAXV];
	  int _vtxHitsOfParticles[MAXV][MAXV];
	  int _tpcHitsOfParticles[MAXV][MAXV];
	  int _ftdHitsOfParticles[MAXV][MAXV];
	  int _sitHitsOfParticles[MAXV][MAXV];
	  int _isProngOfParticles[MAXV][MAXV];
	  float _dEdxOfParticles[MAXV][MAXV];
	  float _errordEdxOfParticles[MAXV][MAXV];
	  int _trueTypeOfParticles[MAXV][MAXV];
	  //TLorentzVector * _particles[MAXV][MAXV]
	  Vertex * _primary;
	  LCCollection * myTrackRel;
	  std::vector<EVENT::MCParticle * > myGenProngs;
	  
	  int _numberOfJets;
	  float _btags[MAXV];
	  float _ctags[MAXV];
	  int _mcpdg[MAXV];
	  int _nvertices[MAXV];
	  float _costhetaJetParticles[MAXV][MAXV3];
	  int _nJetParticles[MAXV];
	  float _vtxAngleJetAxis[MAXV];
	  float _tagAngleJetAxis[MAXV];
	  float _pJetParticles[MAXV][MAXV3];
	  float _zJetParticles[MAXV][MAXV3];
	  float _alphaJetParticles[MAXV][MAXV3];
	  float _maxalphaJetParticles[MAXV];
	  int _typeJetParticles[MAXV][MAXV3];
	  int _prongJetParticles[MAXV][MAXV3];
	  int _vtxJetParticles[MAXV][MAXV3];
	  int _numberOfMissed;
	  float _offsetMissed[MAXV];
	  float _realoffsetMissed[MAXV];
	  int _interactionMissed[MAXV];
	  float _momentumMissed[MAXV];
	  float _massMissed[MAXV];
	  float _thetaMissed[MAXV];
	  float _phiMissed[MAXV];
	  float _ptMissed[MAXV];
	  float _deltapMissed[MAXV];
	  int _genMissed[MAXV];
	  int _interactedMissed[MAXV];
	  float _sigmaAngleMissed[MAXV];
	  float _costhetaMissed[MAXV];
	  int _isrecoMissed[MAXV];
	  int _hastrackMissed[MAXV];
	  float _omegatrackMissed[MAXV];
	  float _distanceIPMissed[MAXV];
	  float _chargeMissed[MAXV];
	  float _vertexAngleMissed[MAXV];
	  float _deviationMissed[MAXV];
	  float _d0Missed[MAXV];
	  float _z0Missed[MAXV];
	  
	  float _truthAngleMissed[MAXV];
	  float _chi2Missed[MAXV];
	  float _errorMissed[MAXV];
	  float _btagMissed[MAXV];
	  int _vtxHitsMissed[MAXV];
	  int _ftdHitsMissed[MAXV];
	  int _tpcHitsMissed[MAXV];
	  int _sitHitsMissed[MAXV];
	
	  int _numberOfMissedVtx;
	  float _costhetaMissedVtx[MAXV];
	  int _genMissedVtx[MAXV];
	  float _momentumMissedVtx[MAXV];
	  float _distanceMissedVtx[MAXV];
	  int _numberOfTracksMissedVtx[MAXV];
	  float _momentumOfParticlesVtx[MAXV][MAXV];
	  float _costhetaOfParticlesVtx[MAXV][MAXV];
	  float _ptOfParticlesVtx[MAXV][MAXV];
	  float _offsetOfParticlesVtx[MAXV][MAXV];

	  int _numberOfSecondaries;
	  int _misrecoNumber;
	  int _primaryNumber;
	  int _wovertexNumber;
	  int _correctNumber;

	  int _numberIncoming;
	  int _numberTagged;
	  int _numberMissed1;
	  int _numberMissed2;
	  int _numberNotHandled;
	} ;
		
} /* TTbarAnalysis */

#endif



