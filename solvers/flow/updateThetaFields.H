/**********************************************/
//-- Update fields related to
//-- unsaturated properties, e.g.
//-- Saturation, Capillary capacity,
//-- Relative permeability, Mobility, fluxes
/**********************************************/
//---------------------------------//
//-- Update saturation
//---------------------------------//
Se_tmp =  max(min(pow(1.0+pow(-min(hp*alpha_vg,-1e-10),n_vg),-m_vg),1-1e-18),1e-18);//for picard hp -> hp_tmp
sw     = sw_min + (sw_max-sw_min)*Se_tmp;

//---------------------------------//
//-- Update Capillary capacity (stor=Ch)
//---------------------------------//
stor.dimensions().reset(dimless/dimLength);
stor = alpha_vg*eps*(sw_max-sw_min)*m_vg/(1.-m_vg)*pow(Se_tmp,1./m_vg)*pow(1.-pow(Se_tmp,1./m_vg),m_vg) ;//*eps needed

//---------------------------------//
//-- Update relative permeability
//---------------------------------//
krw  = pow(Se_tmp,0.5) * pow(1-pow(1-pow(Se_tmp,1/m_vg),m_vg),2); //krmax = 1
krwf = fvc::interpolate(krw,"krw");

//---------------------------------//
//-- Update mobility
//---------------------------------//
Mf = krwf*Kf*rhowf*mag(g)/muwf;	
Lf = krwf*Kf*rhowf/muwf;

//---------------------------------//
//-- Update fluxes
//---------------------------------//
phiGr = (Lf * g) & mesh.Sf();