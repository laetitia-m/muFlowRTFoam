/**********************************************/
//-- Solve the saturation-pressure equation
//-- Variable p and sw
//-- For two phase flow
/**********************************************/

Info << "\n\n//*******************************************************************//" << endl;
Info << "//* SOLVE PRESSURE and SATURATION EQUATIONS (p-Sw) for two-phase flow *//" << endl;
Info << "//*******************************************************************//" << endl;

//---------------------------------//
//-- Variables and initialisation
//---------------------------------//
//-- Time
float time = mesh.time().value();

//-- Source terms
pSu *= 0;
fSu *= 0;

//-- Resize ??
bcPS.resize(nxyz, 0);

//********************************************//
//** Saturation equation 
//********************************************//
Info << "--->> Start solving saturation equation (Sw)..." << endl;
//-- OA: nb as S is already used for solution, we use F instead (fraction of water)

//---------------------------------//
//-- Boundary conditions
//---------------------------------//
Info << "   - Set up the boundary conditions... " << endl;

//-- Saturation well boundary condition
if (cellsFwel.size() > 0)
{
	Info << "     -- Saturation:: well boundary condition" << endl;
	if (time > tFwel)
	{
		output = getCbuffer("fwel", itFwel, cellsFwel.size());
		tFwel = output.t;
		cellsFwelData = output.d;
		itFwel += 1;
		tnext = min(tnext, tFwel);
	}
	for (size_t i = 0; i < cellsFwel.size(); i++)
	{
		i1 = cellsFwel[i];
		qFwel[i] = cellsFwelData[i * 4 + 2];
		fSu[i1] = qFwel[i] / mesh.V()[i1];
		bcPS[i1] = 1;
	}
}

//---------------------------------//
//-- Solve the saturation equation sw
//---------------------------------//
Info << "   - Solve the saturation equation (Sw)... " << endl;
//-- Divergence of water fluxes
volScalarField dphiw = fvc::div(phiw); // fluxes

//-- Saturation equation
fvScalarMatrix SEqn
(
	eps *fvm::ddt(sw) 
	+ fvc::div(phiw) 
	== fSu
);

//-- Get the fixed saturation conditions
// LM:: not clear why this is needed
if (cellsFfix.size() > 0)
{
	if (time > tFfix)
	{
		output = getCbuffer("ffix", itFfix, cellsFfix.size());
		tFfix = output.t;
		cellsFfixData = output.d;
		itFfix += 1;
		tnext = min(tnext, tFfix);
	}
	for (size_t i = 0; i < cellsFfix.size(); i++)
	{
		SEqn.setReference(cellsFfix[i], cellsFfixData[i * 4 + 2]);
		bcPS[cellsFfix[i]] = 1;
	}
}

//-- Solve 
SEqn.solve();

//---------------------------------//
//-- Update of fields (update of what??)
//---------------------------------//
if (plugPS.flag)
{
	#include "plugins/plugin_PSupdate.H"
}
else
{
	#include "flow/update2phaseFields.H"
}

Info << "--->> End solving saturation equation (Sw). " << endl;

//********************************************//
//** Pressure equation 
//********************************************//
Info << "   - Solve the pressure equation (p)... " << endl;

//---------------------------------//
//-- Boundary conditions
//---------------------------------//
Info << "   - Set up the boundary conditions... " << endl;

//-- Pressure well boundary condition
if (cellsPwel.size() > 0)
{
	Info << "     -- Pressure:: well boundary condition" << endl;
	if (time > tPwel)
	{
		output = getCbuffer("pwel", itPwel, cellsPwel.size());
		tPwel = output.t;
		cellsPwelData = output.d;
		itPwel += 1;
		tnext = min(tnext, tPwel);
	}
	for (size_t i = 0; i < cellsPwel.size(); i++)
	{
		i1 = cellsPwel[i];
		qPwel[i] = cellsPwelData[i * 4 + 2];
		pSu[i1] = qPwel[i] / mesh.V()[i1];
		bcPS[i1] = 1;
	}
}

//-- If reaction and boiling
if ((reactIndicator == 1) && (activateEbullition == 1))
{
	Info << "     -- Pressure:: add source term if reaction and boiling" << endl;
	for (j = 0; j < nxyz; j++)
	{
		pSu[j] += (freak.p[j] / (p.oldTime()[j] / atmPa) - 1) * eps[j] * (sw_max[j] - sw.oldTime()[j]) / float(runTime.deltaTValue());
	}
}

//-- Comp_g ??
// comp_g defined as the gas volume gvol divided by the pressure p
// comp_g is never used in the code or the equations
for (j = 0; j < nxyz; j++)
{
	comp_g[j] = gvol[j] / p[j];
}

//---------------------------------//
//-- Solve the pressure equation p
//---------------------------------//
Info << "   - Solve the pressure equation (p)... " << endl;

//-- Pressure equation
fvScalarMatrix pEqn
(
	// comp_g*fvm::ddt(p)
	-fvm::laplacian(Mf, p) 
	+ fvc::div(phiGr) 
	+ fvc::div(phiPc) * activateCapillarity 
	== pSu
);

//-- Get the fixed pressure conditions
// LM:: not clear why this is needed
if (cellsPfix.size() > 0)
{
	if (time > tPfix)
	{
		output = getCbuffer("pfix", itPfix, cellsPfix.size());
		tPfix = output.t;
		cellsPfixData = output.d;
		itPfix += 1;
		tnext = min(tnext, tPfix);
	}
	for (size_t i = 0; i < cellsPfix.size(); i++)
	{
		pEqn.setReference(cellsPfix[i], cellsPfixData[i * 4 + 2]);
		bcPS[cellsPfix[i]] = 1;
	}
}

//-- Solve
pEqn.solve(); // solve

//---------------------------------//
//-- Compute fluxes
//---------------------------------//
Info << "   - Compute and update fluxes... " << endl;

// OA: == to transfer also BC 
// OA: in phiw - keep Fwf (not Mwf/Mf, don't know why)
// OA: in phiw - removed (1-Fwf)*

phiP == -(Mf *fvc::snGrad(p)) * mesh.magSf();
phiw == (Lwf / Lf) * phiGr + Fwf * phiP + phiPc * activateCapillarity; 
phi  == phiGr + phiP + phiPc *activateCapillarity;
phig == phi - phiw;

volScalarField q1 = fvc::div(phiw);

//-- Recompute fluxes at boundaries for fixed BCs
for (size_t i = 0; i < cellsPfix.size(); i++)
{
	i1 = cellsPfix[i];
	qPfix[i] = q1[i1] * mesh.V()[i1];
}

//---------------------------------//
//-- Compute velocities
//---------------------------------//
Info << "   - Compute velocities U, Uw and Ug... " << endl;

U  = fvc::reconstruct(phi);
Uw = fvc::reconstruct(phiw);
Ug = U - Uw;

U.correctBoundaryConditions();
Ug.correctBoundaryConditions();
Uw.correctBoundaryConditions();

Info << "--->> End solving pressure equation (p). " << endl;
Info << "//******************************************************************//\n\n" << endl;
