#include "TrashRecoProcessor.hh"
using std::vector;
using EVENT::Vertex;
using TTbarAnalysis::MathOperator;

namespace TTbarAnalysis 
{
	TrashRecoProcessor aTrashRecoProcessor ;
	TrashRecoProcessor::TrashRecoProcessor() : Processor("TrashRecoProcessor") {
	
	    // modify processor description
	    _description = "TrashRecoProcessor does whatever it does ..." ;
	
	
	    // register steering parameters: name, description, class-variable, default value
	    registerInputCollection( LCIO::VERTEX,
	            "SecondaryCollectionName" , 
	            "Name of the BuildUpVertex collection"  ,
	            _colSecName ,
	            std::string("BuildUpVertex")
	    );
	    registerInputCollection( LCIO::VERTEX,
	            "PrimaryCollectionName" , 
	            "Name of the PrimaryVertex collection"  ,
	            _colPriName ,
	            std::string("PrimaryVertex")
	
	    );
	    registerInputCollection( LCIO::VERTEX,
	    	"MCVertexCollectionName",
		"Name of the MCVertex collection",
		_colMCName,
	    	std::string("MCVertex")
	    );
	}
	
	
	
	void TrashRecoProcessor::init() 
	{ 
		streamlog_out(DEBUG) << "   init called  " << std::endl ;
		printParameters() ;
		_nRun = 0 ;
		_nEvt = 0 ;
		_hfilename = "TrashRecoTest.root";
		_hfile = new TFile( _hfilename.c_str(), "RECREATE", _hfilename.c_str() ) ;
		_hTree = new TTree( "Stats", "My vertex tree!" );
		_hTree->Branch("numberOfTagged", &_numberOfTagged, "numberOfTagged/I");
		_hTree->Branch("numberOfTotal", &_numberOfTotal, "numberOfTotal/I");
		_hTree->Branch("numberOfTernary", &_numberOfTernary, "numberOfTernary/I");
		_hTree->Branch("numberOfSecondary", &_numberOfSecondary, "numberOfSecondary/I");
		_hTree->Branch("numberOfUnknown", &_numberOfUnknown, "numberOfUnknown/I");
		_hTree->Branch("numberOfDistances", &_numberOfDistances, "numberOfDistances/I");
		_hTree->Branch("distances", _distances, "distances[numberOfDistances]/F");
		_hTaggedTree = new TTree( "TaggedVertices", "My vertex tree!" );
		_hTaggedTree->Branch("numberOfTagged", &_numberOfTagged, "numberOfTagged/I");
		_hTaggedTree->Branch("distance", _distanceFromIP, "distance[numberOfTagged]/F");
		_hTaggedTree->Branch("mindistance", _mindistance, "mindistance[numberOfTagged]/F");
		_hTaggedTree->Branch("probability", _probability, "probability[numberOfTagged]/F");
		_hTaggedTree->Branch("chi2", _chi2, "chi2[numberOfTagged]/F");
		_hTaggedTree->Branch("coordinates", _coordinates, "coordinates[numberOfTagged][3]/F");
		_hTaggedTree->Branch("PDG", _PDG, "PDG[numberOfTagged]/I");
		_hTaggedTree->Branch("generation", _generation, "generation[numberOfTagged]/I");
		_hTaggedTree->Branch("charge", _charge, "charge[numberOfTagged]/I");
		_hTaggedTree->Branch("numberOfParticles", _numberOfParticles, "numberOfParticles[numberOfTagged]/I");
		_hTaggedTree->Branch("energyOfParticles", _energyOfParticles, "energyOfParticles[numberOfTagged][15]/F");
		_hTaggedTree->Branch("momentumOfParticles", _momentumOfParticles, "momentumOfParticles[numberOfTagged][15]/F");
		_hTaggedTree->Branch("massOfParticles", _massOfParticles, "massOfParticles[numberOfTagged][15]/F");
		
		_hUntaggedTree = new TTree( "UntaggedVertices", "My vertex tree!" );
		_hUntaggedTree->Branch("numberOfUnknown", &_numberOfUnknown, "numberOfUnknown/I");
		_hUntaggedTree->Branch("distance", _distanceFromIP, "distance[numberOfUnknown]/F");
		_hUntaggedTree->Branch("probability", _probability, "probability[numberOfUnknown]/F");
		_hUntaggedTree->Branch("chi2", _chi2, "chi2[numberOfUnknown]/F");
		_hUntaggedTree->Branch("charge", _charge, "charge[numberOfUnknown]/I");
		_hUntaggedTree->Branch("coordinates", _coordinates, "coordinates[numberOfUnknown][3]/F");
		_hUntaggedTree->Branch("numberOfParticles", _numberOfParticles, "numberOfParticles[numberOfUnknown]/I");
		_hUntaggedTree->Branch("energyOfParticles", _energyOfParticles, "energyOfParticles[numberOfUnknown][15]/F");
		_hUntaggedTree->Branch("momentumOfParticles", _momentumOfParticles, "momentumOfParticles[numberOfUnknown][15]/F");
		_hUntaggedTree->Branch("massOfParticles", _massOfParticles, "massOfParticles[numberOfUnknown][15]/F");
		//_hTaggedTree->Branch("particles", _particles);
		
	
	}
	
	
	void TrashRecoProcessor::processRunHeader( LCRunHeader* run) 
	{ 
		_nRun++ ;
	} 
	
	Vertex * TrashRecoProcessor::FindPrimaryVertex(const LCCollection * col)
	{
		int number = col->getNumberOfElements();
		if (number > 1) 
		{
			std::cout << "WARNING: Pile-up found!\n";
		}
		return dynamic_cast< Vertex * >( col->getElementAt(0) ) ;
	}
	void TrashRecoProcessor::PrintParticle(ReconstructedParticle * particle)
	{
		if (!particle) 
		{
			return;
		}
		std::cout << std::fixed << std::setw( 6 ) << std::setprecision( 3 ) << std::setfill( ' ' );
		int id = 0;
		if (particle->getParticleIDUsed()) 
		{
			std::cout << "Type " << particle->getParticleIDUsed()->getType() << '\n';
			id = particle->getParticleIDs()[0]->getPDG(); 
		}
		std::cout<<"|"<< id <<"\t\t|"<<particle->getMass()<<"\t\t|"<<particle->getCharge()  <<"\t\t|"<<particle->getEnergy() <<"\t\t|\n";
	
	}
	
	void TrashRecoProcessor::processEvent( LCEvent * evt ) 
	{ 
		try
		{
			LCCollection* col = evt->getCollection( _colSecName );
			LCCollection* mc = evt->getCollection( _colMCName );
			int number = col->getNumberOfElements();
			VertexRecoOperator reco;
			//vector< VertexTag * > * tagged = reco.Compare(col, mc);
			vector< VertexTag * > * tagged = reco.CompareDirection(col, mc);
			_primary = FindPrimaryVertex( evt->getCollection( _colPriName ));
			if (_primary) 
			{
				_numberOfTagged = tagged->size();
				_numberOfUnknown = number - _numberOfTagged;
				std::cout << "We have " << number << " vertices\n";
				for (int i = 0; i < _numberOfTagged; i++) 
				{
					Vertex * vertex = tagged->at(i)->GetVertex();
					ReconstructedParticle * particle = vertex->getAssociatedParticle();
					PrintParticle(particle);
					for (int j = 0; j < particle->getParticles().size(); j++) 
					{
						PrintParticle(particle->getParticles()[j]);
					}
					Write(tagged->at(i), i);
				}
				for (int i = 0; i < number; i++) 
				{
					Vertex * vertex = dynamic_cast< Vertex * >( col->getElementAt(i) ) ;
				}
				
			}
			_numberOfTotal = number;
			/*vector<float> distances = reco.GetDistances();
			_numberOfDistances = distances.size();
			std::cout << "N_dist: " << _numberOfDistances << '\n';
			for (int i = 0; i < _numberOfDistances; i++) 
			{
				_distances[i] = distances[i];
			}*/
			_hTree->Fill();
			_hTaggedTree->Fill();
			ClearVariables();
			vector< Vertex * > unknown = reco.GetUnknownVertexes();
			for (int i = 0; i < unknown.size(); i++) 
			{
				Write(unknown[i],i);
			}
			std::cout << "We have " << unknown.size() << " unknown vertices\n";
			_numberOfUnknown = unknown.size();
			_hUntaggedTree->Fill();
			ClearVariables();
	
		}
		catch( DataNotAvailableException &e)
		{
			streamlog_out(DEBUG) << "No collection!" << std::endl ;
		}
		_nEvt ++ ;
	}
	void TrashRecoProcessor::Write (Vertex * vertex, int number)
	{
		const float * position = vertex->getPosition();
		for (int i = 0; i < 3; i++) 
		{
			_coordinates[number][i] = position[i];
		}
		_probability[number] = vertex->getProbability();
		_chi2[number] = vertex->getChi2();
		ReconstructedParticle * particle = vertex->getAssociatedParticle();
		_charge[number] = particle->getCharge();
		_numberOfParticles[number] = particle->getParticles().size();
		for (int j = 0; j < _numberOfParticles[number]; j++)
		{
			ReconstructedParticle * component = particle->getParticles()[j];
			_momentumOfParticles[number][j] = MathOperator::getModule( component->getMomentum());
			_energyOfParticles[number][j] = component->getEnergy();
			_massOfParticles[number][j] = component->getMass();
		}
		float distance = MathOperator::getDistance(_primary->getPosition(), vertex->getPosition() );
		_distanceFromIP[number] = distance;
		std::cout << "d(V_p,V_i) = " << distance << "; P(V_i) = " << vertex->getProbability() << "; Chi^2(V_i) = " << vertex->getChi2() << '\n';
		
	}
	void TrashRecoProcessor::Write (VertexTag * tag, int number)
	{
		if (!tag) 
		{
			return;
		}
		_numberOfTernary += (tag->GetGeneration() == 3)? 1 : 0;
		_numberOfSecondary += (tag->GetGeneration() == 2)? 1 : 0;
		_PDG[number] = tag->GetInitialPDG();
		_generation[number] = tag->GetGeneration();
		_mindistance[number] = tag->GetMinimalDistance();
		Vertex * vertex = tag->GetVertex();
		Write (vertex, number);
	}
	void TrashRecoProcessor::ClearVariables()
	{
		_numberOfTernary = 0;
		_numberOfSecondary = 0;
		_numberOfTagged = 0;
		_numberOfTotal = 0;
		_numberOfUnknown = 0;
		for (int i = 0; i < MAXV2; i++) 
		{
			_distances[i] = -1.0;
		}
		for (int i = 0; i < MAXV; i++) 
		{
			_probability[i] = -1.0;
			_PDG[i] = 0;
			_generation[i] = 0;
			_chi2[i] = -1.0;
			_distanceFromIP[i] = -1.0;
			for (int j = 0; j < MAXV; j++) 
			{
				_energyOfParticles[i][j] = -1.0;
				_momentumOfParticles[i][j] = -1.0;
				_massOfParticles[i][j] = -1.0;
			}
		}
	}
	
	void TrashRecoProcessor::check( LCEvent * evt ) 
	{ 
		// nothing to check here - could be used to fill checkplots in reconstruction processor
		
	}
	
	
	void TrashRecoProcessor::end()
	{ 
		_hfile->cd();
		_hfile->Write();
		_hfile->Close();
	    //   std::cout << "TrashRecoProcessor::end()  " << name() 
	    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    // 	    << std::endl ;
	
	}
} /* TTbarAnalysis */
