
/**********************************************/
//-- Update fields related to
//-- two phases, e.g.
//-- Saturation, Capillary pressure,
//-- Permeabilities, Mobilities, fluxes
/**********************************************/
/**********************************************/
double mysmall = 1e-18;

//---------------------------------//
//-- Update saturation
//---------------------------------//
sw = sw_max - sg;
Se = min(max((sw-sw_min)/(sw_max-sw_min),mysmall),1-mysmall); // Effective water saturation

//---------------------------------//
//-- Update capillary pressure and gradient of Pc
//---------------------------------//
pc == pc0*pow(pow(Se,-1/m_vg)-1,1/n_vg); // == instead of =? Should return a boolean?
gradpc  = fvc::grad(pc,"pc");//
gradpcf = fvc::interpolate(gradpc,"gradpc");//

//---------------------------------//
//-- Update permeabilities
//---------------------------------//
krw  = 1 * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_vg),m_vg),2);
krg  = pow(1-Se,0.5) * pow(1-pow(Se,1/m_vg),2*m_vg);
krg  = max(min(krg,1-mysmall),mysmall); 

krwf = fvc::interpolate(krw,"krw");
krgf = fvc::interpolate(krg,"krg");

Fwf  = (krwf/muwf) / (krgf/mug + krwf/muwf);

//---------------------------------//
//-- Update mobilities and fractional flow
//---------------------------------//
Mwf = Kf *krwf/muwf;
Mgf = Kf *krgf/mug;
Mf  = Mwf + Mgf;

Lwf = Mwf*rhowf;
Lgf = Mgf*rhog;
Lf  = Lwf + Lgf;

//---------------------------------//
//-- Update Fluxes
//---------------------------------//
phiGr = (Lf * g) & mesh.Sf();
phiPc = (Mwf * gradpcf) & mesh.Sf(); 