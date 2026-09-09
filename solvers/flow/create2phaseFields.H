/**********************************************/
//-- Creation of fields related to
//-- unsaturated properties
/**********************************************/
Info << "Reading and creating 2 phase fields and properties... FlowType is " << flowType << endl;

//---------------------------------//
//-- Density and viscosity of gas
//---------------------------------//
dimensionedScalar rhog
(
	transportProperties.lookupOrDefault
	(
		"rhog",
		dimensionedScalar("rhog",dimensionSet(1,-3,0,0,0,0,0),1.2) // Reference value of air
	)
);

dimensionedScalar mug
(
	transportProperties.lookupOrDefault
	(
		"mug",
		dimensionedScalar("mug",dimensionSet(1,-1,-1,0,0,0,0),1e-5) // Reference value of air
	)
);

//---------------------------------//
//-- Gas Velocity
//---------------------------------//
volVectorField Ug
(
    IOobject
	(
		"Ug",
		TIME_NAME(mesh.time()),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedVector("Ug",dimLength/dimTime,vector::zero)
);

//-- Total velocity
volVectorField U = Ug + Uw;

// surfaceScalarField phig ("phig",phiw*0); // LM: I do not understand why this line needs to be exactly here

//---------------------------------//
//-- Pressure
//---------------------------------//
volScalarField p1=p;
volScalarField pw=p;
volScalarField p_tmp=p;

//---------------------------------//
//-- Gas Saturation
//---------------------------------//
volScalarField sg 
(
	IOobject
	(
		"sg",
		TIME_NAME(mesh.time()),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("sg",dimless,0)
);

//-- Effective gas saturation
Se = max(SMALL,min((sw-sw_min)/(sw_max-sw_min),1.-SMALL)); // Effective water saturation
volScalarField Seg = Se; // Effective gas saturation

//---------------------------------//
//-- comp_g ??
//---------------------------------//
// comp_g is defined as the gas volume gvol divided by the pressure p
// comp_g=gvol/p
// However, comp_g is never used in the code or the equations
volScalarField comp_g
(
	IOobject
	(
		"comp_g",
		runTime.constant(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("comp_g",dimensionSet(-1,1,2,0,0,0,0),1./atmPa)
);

//---------------------------------//
//-- Fluxes 
//---------------------------------//
//-- Fluxes are computed here, because it is needed in interpolation function
surfaceScalarField phi = phiw;
surfaceScalarField phig ("phig",phiw*0); // LM: This flux needs to be here because it is used in fvc:interpolate
surfaceScalarField phiPc ("phiPc",phiw*0); // LM: use to have an equal phiPc = ("phiPc",phiw*0), not sure what is the difference...
surfaceScalarField phiP ("phiP", phiw);

//---------------------------------//
//-- Relative permeabilities of gas
//---------------------------------//
volScalarField krg ("krg", pow(1-Se,0.5) * pow(1-pow(Se,1/m_vg),2*m_vg) );
surfaceScalarField krgf ("krgf",fvc::interpolate(krg,"krg"));

//---------------------------------//
//-- Capillary pressure and gradient of Pc
//---------------------------------//
volScalarField pc
(
	IOobject
	(
		"pc",
		runTime.constant(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("pc",dimensionSet(1,-1,-2,0,0,0,0),1.)
);

volScalarField pc0 = mag(g)*rhow/alpha_vg; // be careful, if g=0 pc0=0
pc==pc0*pow(pow(Se+VSMALL,-1./m_vg)-1,1./n_vg); // == instead of =? Should return a boolean?

// Gradient of capillary pressure
volVectorField gradpc ("gradpc", fvc::grad(pc,"pc") );
surfaceVectorField gradpcf ("gradpcf", fvc::interpolate(gradpc,"gradpc") );

// Derivative of capillary pressure with respect to water saturation
volScalarField dpcdS ("dpcdS", pc0*0);
surfaceScalarField dpcdSf ("dpcdSf", fvc::interpolate(dpcdS) );

//---------------------------------//
//-- Mobility Mf and Fractional flow Lf
//---------------------------------//
surfaceScalarField Mwf ("Mwf",Kf*krwf/muwf);
surfaceScalarField Mgf ("Mgf",Kf*krgf/mug);

surfaceScalarField Lwf ("Lwf",Mwf*rhowf);
surfaceScalarField Lgf ("Lgf",Mgf*rhog);

surfaceScalarField Fwf ("Fwf",(krwf/muwf)/(krgf/mug + krwf/muwf));

if (flowType == 4) // 2 phase
{
	Mf.dimensions().reset(dimensionSet(-1,3,1,0,0,0,0));
	Lf.dimensions().reset(dimensionSet(0,0,1,0,0,0,0));

	Mf = Mwf + Kf*krgf/mug;
	Lf = Lwf + Kf*rhog*krgf/mug;
}

//---------------------------------//
//-- Fluxes for 2-phase
//---------------------------------//
if (flowType == 4) // 2 phase
{
	surfaceScalarField phiPc ("phiPc", (Mwf * gradpcf) & mesh.Sf() );
	surfaceScalarField phiGr ("phiGr", (Lf * g) & mesh.Sf());
}

//---------------------------------//
//-- b ?? Seems to not be used anymore
//---------------------------------//
if (flowType == 4) // 2 phase
{
	b.dimensions().reset(dimensionSet(-1,1,2,0,0,0,0));
}

//---------------------------------//
//-- Source terms - pSu and fSu
//---------------------------------//
// Used in pSEqn.H

// pSu: Source term of p equation
// pSu = qPwel/cell_vol
volScalarField::Internal pSu
(
   IOobject
   (
		"pSu",
		TIME_NAME(mesh.time()),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0) // the overall eq is in m3.s-1
);

// fSu: Source term of Sw equation
// fSu = qFwel/cell_vol
volScalarField::Internal fSu
(
   IOobject(
		"fSu",
		TIME_NAME(mesh.time()),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0),0.0) // the overall eq is in m3.s-1
);

//---------------------------------//
//-- Condition for pressure - LM: whatever this means
//---------------------------------//
//-- Declare variables
labelList cellsPfix0; std::vector<float> qPfix,cellsPfixData;
labelList cellsFfix0; std::vector<float> qFfix,cellsFfixData; // LM: not used
labelList cellsPwel0; std::vector<float> qPwel,cellsPwelData;
labelList cellsFwel0; std::vector<float> qFwel,cellsFwelData;
int   itPfix, itFfix, itPwel, itFwel;
float tPfix,  tFfix,  tPwel,  tFwel;
std::vector<int>  cellsPfix,cellsPwel,cellsFfix,cellsFwel;

// Not sure what this is for. Seems to assign 1 to cells with one the conditions above
std::vector<int> bcPS(ncell,0);

//-- For qPfix
if (fexists(cur_dir+"/constant/options/pfix"))
{ 
	std::cout<<"in pfix reader \n";
	cellSet selectedPfix(mesh, "pfix");
	cellsPfix0 = selectedPfix.toc();
	qPfix.resize(cellsPfix0.size());
	while (time>=tPfix) 
	{
		output = getCbuffer("pfix",itPfix,cellsPfix0.size());
		itPfix +=1;
		tPfix = output.t;
		cellsPfixData=output.d;
	}
	cellsPfix = indexC(cellsPfix0,cellsPfixData) ;
}

//-- For qFfix - LM: but this seems to not be used anywhere.
if (fexists(cur_dir+"/constant/options/ffix")) 
{ 
	cellSet selectedFfix(mesh, "ffix");
	cellsFfix0 = selectedFfix.toc();
	qFfix.resize(cellsFfix0.size());
	while (time>=tFfix) 
	{
		output = getCbuffer("ffix",itFfix,cellsFfix0.size());
		itFfix +=1;
		tFfix = output.t;
		cellsFfixData=output.d;
	}
	cellsFfix = indexC(cellsFfix0,cellsFfixData) ;
}

//-- For qPwel
if (fexists(cur_dir+"/constant/options/pwel")) 
{ 
	cellSet selectedPwel(mesh, "pwel");
	cellsPwel0 = selectedPwel.toc();
	qPwel.resize(cellsPwel0.size());
	while (time>=tPwel) 
	{
		output = getCbuffer("pwel",itPwel,cellsPwel0.size());
		itPwel +=1;
		tPwel = output.t;
		cellsPwelData=output.d;
	}
	cellsPwel = indexC(cellsPwel0,cellsPwelData) ;
}

//-- For qFwel
if (fexists(cur_dir+"/constant/options/fwel")) 
{ 
	cellSet selectedFwel(mesh, "fwel");
	cellsFwel0 = selectedFwel.toc();
	qFwel.resize(cellsFwel0.size());
	while (time>=tFwel) 
	{
		output = getCbuffer("fwel",itFwel,cellsFwel0.size());
		itFwel +=1;
		tFwel = output.t;
		cellsFwelData=output.d;
	}
	cellsFwel = indexC(cellsFwel0,cellsFwelData) ;
}
