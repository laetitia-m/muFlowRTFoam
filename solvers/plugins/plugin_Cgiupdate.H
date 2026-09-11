//change the values used for Cgi equation
for (size_t i=0;i<listCouples.size();i++) {
	volScalarField X = mesh.objectRegistry::lookupObject<volScalarField>(listCouples[i].Xvar);  //Xvar is the name of the variable

	//--->> START:: Added by ZA
	volScalarField X2 = mesh.objectRegistry::lookupObject<volScalarField>(listCouples[i].Xvar2);
	volScalarField X3 = mesh.objectRegistry::lookupObject<volScalarField>(listCouples[i].Xvar3);
	//<<-- Added by ZA:: END
	
	Info<<" X name "<<listCouples[i].Xvar<<" X(0) "<<X[0]<<endl;
	if (listCouples[i].Ykey=="Deffg")
		{
		//--->> START:: Added by ZA
		// Deffg *= calcFunc(mesh,listCouples[i].type,X,listCouples[i].Xref,listCouples[i].fparms);
		Deffg *= calcFunc(mesh,listCouples[i].type,X,X2,X3,listCouples[i].Xref,listCouples[i].fparms,tunits);
		//<<-- Added by ZA:: END
		Info<<"Deffg "<<Deffg[0]<<endl;
		}
}
