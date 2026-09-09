/*This plugin will modify some elements before the H equation
including transformation of muw, rhow and k
the equations are read from a unique file called coupling and what is sent back is a multiplying factor 
for the given variable
the file coupling is read to find variables of interest (Y) and the predicting variable (X)
this is stored a series of Y=f(X) fnctions that will be applied when needed
*/
class plugin_H
{
public:
	//volScalarField factM; don't put that here, it creates pb 
	std::vector<dicFunc> listCouples;

	//--->> START:: Added by ZA
	std::vector<float> AntonioParms;
	scalar muwRefValue;
	//<<-- Added by ZA:: END

	float tunits;
	
	void init(std::string & cur_dir,IOdictionary & Tprop,fvMesh & mesh,Time & runTime,std::vector<dicFunc> & listCouples,float tunits)
	{
		std::cout<<"init plugin H "<<" ";
		this->listCouples = listCouples;

		//--->> START:: Added by ZA
		dimensionedScalar muwRef(Tprop.lookupOrDefault("muwRef",dimensionedScalar("muwRef",dimensionSet(1,-1,-1,0,0,0,0),1e-3)));
		this->muwRefValue = muwRef.value();
		if (Tprop.found("AntonioCoeffs"))
		{
			const dictionary& coeffs = Tprop.subDict("AntonioCoeffs");
			this->AntonioParms.resize(6);
			this->AntonioParms[0] = readScalar(coeffs.lookup("k0"));
			this->AntonioParms[1] = readScalar(coeffs.lookup("k1"));
			this->AntonioParms[2] = readScalar(coeffs.lookup("dn"));
			this->AntonioParms[3] = readScalar(coeffs.lookup("cN"));
			this->AntonioParms[4] = readScalar(coeffs.lookup("mu1"));
			this->AntonioParms[5] = readScalar(coeffs.lookup("mu2"));
		}
		for (size_t i=0;i<this->listCouples.size();i++)
		{
			if ((this->listCouples[i].type=="Antonio") && (this->listCouples[i].fparms.size()==0))
			{
				if (this->AntonioParms.size() != 6)
				{
					FatalErrorIn("plugin_H::init")
						<< "Antonio selected in coupling with 0 inline parameters, "
						<< "but transportProperties has no complete AntonioCoeffs dictionary"
						<< exit(FatalError);
				}
				this->listCouples[i].setParms(this->AntonioParms);
			}
		}
		//<<-- Added by ZA:: END
	
		this->tunits = tunits;
		//reps0koz = pow(eps,3)/pow(1-eps,2);
	}
	
	volScalarField corrM(fvMesh & mesh,my_phq & freak) //this may modify muw, rhow due to effects of C or T, and K due to minerals (Mf=rhow*mag(g)*Kf/muw)
	{
		volScalarField fact=mesh.objectRegistry::lookupObject<volScalarField>("unity");
		volVectorField Uw =mesh.objectRegistry::lookupObject<volVectorField>("Uw");
		//dimensionedScalar tunits = mesh.objectRegistry::lookupObject<dimensionedScalar>("tunits");
		volScalarField mUw = mag(Uw);
		
		for (size_t i=0;i<listCouples.size();i++) {
			std::cout<<"cpl key "<<listCouples[i].Ykey<<" type "<< listCouples[i].type <<"\n";
			volScalarField X = mesh.objectRegistry::lookupObject<volScalarField>(listCouples[i].Xvar);  //Xvar is the name of the variable

			//--->> START:: Added by ZA
			volScalarField X2 = mesh.objectRegistry::lookupObject<volScalarField>(listCouples[i].Xvar2);  //Xvar is the name of the variable
			volScalarField X3 = mesh.objectRegistry::lookupObject<volScalarField>(listCouples[i].Xvar3);
			//<<-- Added by ZA:: END

			if (listCouples[i].Xvar=="eps") 
				{
				volScalarField eps0 = mesh.objectRegistry::lookupObject<volScalarField>("eps0");  
				volScalarField eps = mesh.objectRegistry::lookupObject<volScalarField>("eps");  
				forAll(eps(),i){eps[i]=eps0[i]-freak.deltaMvol[i]*0.037;} // !! this is valid only for calcite //if (i<5) {Info<<"dm "<<freak.deltaMvol[i]<<endl;}
				X = eps*1;
				} 
			if ((listCouples[i].Ykey=="rhow")||(listCouples[i].Ykey=="K"))
				//--->> START:: Added by ZA
				{fact *= calcFunc(mesh,listCouples[i].type,X,X2, X3,listCouples[i].Xref, listCouples[i].fparms);}
				// {fact *= calcFunc(mesh,listCouples[i].type,X,listCouples[i].Xref,listCouples[i].fparms);}
				//<<-- Added by ZA:: END
				if (listCouples[i].Ykey=="muw")
				//--->> START:: Added by ZA
				{fact /=calcFunc(mesh,listCouples[i].type,X,X2,X3,listCouples[i].Xref, listCouples[i].fparms);}
				// {fact = 1.0/calcFunc(mesh,listCouples[i].type,X,listCouples[i].Xref,listCouples[i].fparms);}
				//<<-- Added by ZA:: END
				if (listCouples[i].Ykey=="muShear") 
				{fact = 1.0/calcFunc2(mesh,listCouples[i].type,X,mUw,listCouples[i].Xref,listCouples[i].fparms,tunits);}
		}
		
		return fact;
	}

};