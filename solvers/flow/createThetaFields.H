/**********************************************/
//-- Creation of fields related to
//-- unsaturated properties
/**********************************************/
Info << "Reading unsaturated properties... FlowType is " << flowType << endl;

//---------------------------------//
//-- Pressure head hp
//---------------------------------//
Info << "Reading pressure head hp.\n" << endl;
volScalarField hp
(
    IOobject
	(
		"hp",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("hp",dimLength,0)
);
volScalarField hp1=hp;
volScalarField hp_tmp=hp;

//---------------------------------//
//-- Gas and water content
//---------------------------------//
volScalarField thetag ("thetag",unity);
volScalarField thetaw ("thetaw",unity);

//---------------------------------//
//-- Van Genuchten parameters
//---------------------------------//
Info << "Reading Van Genuchten parameters (n_vg, alpha_vg). \n" << endl;

volScalarField n_vg
(
    IOobject
	(
		"n_vg",
		runTime.constant(),
		mesh,IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("n_vg",dimless,1.5)
);

volScalarField alpha_vg
(
    IOobject
	(
		"alpha_vg",
		runTime.constant(),
		mesh,IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("alpha_vg",dimless/dimLength,1.)
);

volScalarField m_vg = 1 - 1/n_vg;

//---------------------------------//
//-- Water Saturation
//---------------------------------//
Info << "Reading min and max water saturation (sw_min, sw_max). \n" << endl;

volScalarField sw_min
(
    IOobject
	(
		"sw_min",
		runTime.constant(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("sw_min",dimless,1e-6)
);

volScalarField sw_max
(
    IOobject
	(
		"sw_max",
		runTime.constant(),
		mesh,IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("sw_max",dimless,1-1e-4)
);

//-- Effective water saturation for kr and Pc Van Genuchten model
volScalarField Se ("Se",unity);
volScalarField Se_tmp ("Se_tmp",unity); // Used in updateThetaFields.H

if (flowType == 3) // unsaturated flow type
{
	Se = pow(1.0+pow(-min(hp*alpha_vg,-1e-10),n_vg),-m_vg)-SMALL; // SMALL is defined by OpenFOAM as being 1e-15
	sw = sw_min + (sw_max-sw_min)*Se;
	Info<<"n_vg "<<min(n_vg)<<" alpha "<<min(alpha_vg)<<" min Se "<<min(Se)<<" max Se "<<max(Se)<<endl;
}

//---------------------------------//
//-- Storage coefficient - Water retention capacity
//---------------------------------//
if (flowType == 3) // unsaturated flow type
{
	// HERE Ch is called stor for compatibility
	stor.dimensions().reset(dimless/dimLength);
	stor = eps*alpha_vg*m_vg*(sw_max-sw_min)/(1-m_vg)*pow(Se,1/m_vg)*pow(1-pow(Se,1/m_vg),m_vg);
}

//---------------------------------//
//-- Relative permeabilities
//---------------------------------//
volScalarField krw ("krw",unity);
surfaceScalarField krwf ("krwf",fvc::interpolate(krw,"krw"));

// Derivative of permeability with respect to saturation.
// Used in CoatsNo.H.
volScalarField dkrwdS ("dkrwdS", unity); 

if (flowType == 3) // unsaturated
{
	krw = 1 * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_vg),m_vg),2);
	krwf = fvc::interpolate(krw,"krw");
}

//---------------------------------//
//-- Mobility Mf and Fractional flow Lf
//---------------------------------//
surfaceScalarField Lf = rhowf*Kf*krwf/muwf;

if (flowType == 3) // unsaturated flow type
{
	Mf = rhowf*mag(g)*Kf*krwf/muwf;// Mf alread exists
	Lf = rhowf*Kf*krwf/muwf; // I think we do not need this as it also defined outside the flowType condition.
}

//---------------------------------//
//-- Fluxes
//---------------------------------//
surfaceScalarField phiGr ("phiGr",(Lf * g) & mesh.Sf() );

if (flowType == 3) // unsaturated flow type
{
	phiGr = (Lf * g) & mesh.Sf();
}

Info << "... end of reading unsaturated properties." <<endl;
