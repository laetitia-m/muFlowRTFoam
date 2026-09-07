/**********************************************/
//-- Solve the diffusion-dispersion equation
//-- For concentration of components in water
//-- Cwi with reactive transport
/**********************************************/

Info << "\n\n//******************************************************************//" << endl;
Info << "//* * * * * * SOLVE DIFFUSION-DISPERSION EQUATION Cwi * * * * * * * * //" << endl;
Info << "//******************************************************************//" << endl;
Info << "--->> Start solving diffusion-dispersion equation (reactive transport - concentration Cwi) " << endl;

//---------------------------------//
//-- Variables
//---------------------------------//
int ic, jmax;
size_t i;

dimensionedScalar USMALL("", dimLength / dimTime, 1e-36);

//- DeltaC max
const scalar dCmax = runTime.controlDict().lookupOrDefault("dCmax", 0.02);

//-- Water content
thetaw = eps * sw;

//-- ??
mxCv.resize(ph_ncomp - 4, 0.);
bcCwi.resize(ncell, 0);

//-- Time
time = mesh.time().value(); // nb tnext is the future time that a new condition will appear

//---------------------------------//
//-- Compute the effective water diffusivity
//---------------------------------//
Deffw = (alphaT * mag(Uw) + Dw0 * pow(thetaw, 3.333) / pow(eps, 2.)) * tensor::I + (alphaL - alphaT) * (Uw * Uw) / (mag(Uw) + USMALL);
Info << "   - Diffusivity: " << Dw0.value() << endl;
Info << "   - Effective diffusivity matrix: " << Deffw[15] << endl;

//---------------------------------//
//-- Boundary conditions
//---------------------------------//
Info << "   - Set up the boundary conditions... " << endl;

if (flowType < 4) //-- If flow is no flow, confined, unconfined or unsaturated
{
	//-- No solutions for drains

	//-- Concentration well boundary condition
	if (cellsSghb.size() > 0)
	{
		Info << "     -- Concentration:: GHB boundary condition" << endl;
		if (time > tSghb)
		{
			output = getCbuffer("sghb", itSghb, cellsSghb.size());
			tSghb = output.t;
			cellsSghbData = output.d;
			itSghb += 1;
			tnext = min(tnext, tSghb);
		}
		for (i = 0; i < cellsSghb.size(); i++)
		{
			condSghb[i] = static_cast<int>(cellsSghbData[i * 4 + 2]);
			i1 = cellsSghb[i];
			bcCwi[i1] = 1;
		}
	}

	//-- Concentration fixed head boundary condition
	if (cellsShfix.size() > 0)
	{
		Info << "     -- Concentration:: fixed head boundary condition" << endl;
		if (time > tShfix)
		{
			output = getCbuffer("shfix", itShfix, cellsShfix.size());
			tShfix = output.t;
			cellsShfixData = output.d;
			itShfix += 1;
			tnext = min(tnext, tShfix);
		}
		for (i = 0; i < cellsShfix.size(); i++)
		{
			condShfix[i] = static_cast<int>(cellsShfixData[i * 4 + 2]);
			i1 = cellsShfix[i];
			bcCwi[i1] = 1;
		}
	}

	//-- Concentration recharge boundary condition
	//-- OA: here only solution given in rch zone)
	if (cellsSrch.size() > 0)
	{
		Info << "     -- Concentration:: recharge boundary condition" << endl;
		if (time > tSrch)
		{
			output = getCbuffer("srch", itSrch, cellsSrch.size());
			tSrch = output.t;
			cellsSrchData = output.d;
			itSrch += 1;
			tnext = min(tnext, tSrch);
		}
		for (i = 0; i < cellsSrch.size(); i++)
		{
			condSrch[i] = static_cast<int>(cellsSrchData[i * 4 + 2]);
			i1 = cellsSrch[i];
			bcCwi[i1] = 1;
		}
	}
}
else // If flowType is 4, i.e. dual phase
{
	//-- Concentration fixed pressure boundary condition
	if (cellsSpfix.size() > 0)
	{
		Info << "     -- Concentration:: fixed pressure boundary condition" << endl;
		if (time > tSpfix)
		{
			output = getCbuffer("spfix", itSpfix, cellsSpfix.size());
			tSpfix = output.t;
			cellsSpfixData = output.d;
			itSpfix += 1;
			tnext = min(tnext, tSpfix);
		}
		for (i = 0; i < cellsSpfix.size(); i++)
		{
			condSpfix[i] = static_cast<int>(cellsSpfixData[i * 4 + 2]);
			i1 = cellsSpfix[i];
			bcCwi[i1] = 1;
		}
	}
}

//-- Concentration fixed boundary condition
if (cellsSfix.size() > 0)
{
	Info << "     -- Concentration:: fixed boundary condition" << endl;
	if (time > tSfix)
	{
		output = getCbuffer("sfix", itSfix, cellsSfix.size());
		tSfix = output.t;
		cellsSfixData = output.d;
		itSfix += 1;
		tnext = min(tnext, tSfix);
	}
	for (i = 0; i < cellsSfix.size(); i++)
	{
		condSfix[i] = static_cast<int>(cellsSfixData[i * 4 + 2]);
		i1 = cellsSfix[i];
		bcCwi[i1] = 1;
	}
}

//-- Concentration well boundary condition
if (cellsSwel.size() > 0)
{
	if (time > tSwel)
	{
		Info << "     -- Concentration:: well boundary condition" << endl;
		output = getCbuffer("swel", itSwel, cellsSwel.size());
		tSwel = output.t;
		cellsSwelData = output.d;
		itSwel += 1;
		tnext = min(tnext, tSwel);
	}
	for (i = 0; i < cellsSwel.size(); i++)
	{
		condSwel[i] = static_cast<int>(cellsSwelData[i * 4 + 2]);
		i1 = cellsSwel[i];
		bcCwi[i1] = 1;
	}
}

//---------------------------------//
//-- Loop over species
//---------------------------------//
float dC = 0;
dtForC = 1e12; // Time-step computed for each species

//-- Main loop over species
unity = unity * 0 + 1;
inactiveFaces.clear();

//-- OA: water is always set to 55.34 from solution 0
for (j = 0; j < nxyz; j++)
{
	Cw[0]()[j] = solu_conc[0];
} 

//-- Loop over species
for (ic = 1; ic < ph_ncomp; ic++)
{
	volScalarField &Cwi = Cw[ic]();
	word nameCi("Cw" + name(ic));
	if (immobile[ic] == 1)
	{
		continue;
	}

	//-- Previous cells with values must be destroyed
	CiSu *= 0;
	CiSp *= 0;

	// --------------------- internal conditions ------------------------
	// we add phreeqc results as source term (during one tstep) : CiSu : M/time/Volume, cellsHwelas C=M/V -> Csu=(C_pqh-C_opf)/dt
	// these results will be replaced by classical ones at well and hfix bc (no chem at wells)
	Cw[ic]().storePrevIter();

	//---------------------------------//
	//-- Boundary conditions
	//---------------------------------//
	if (reactIndicator == 1)
	{
		for (j = 0; j < nxyz; j++)
		{
			CiSu[j] += (Cw[ic]()[ractive[j]] - freak.c[ic * nxyz + j]) / float(runTime.deltaTValue()) * thetaw[j];
		}
	}

	if (flowType < 4) //-- If flow is no flow, confined, unconfined or unsaturated
	{
		//-- Head:: well boundary condition
		//-- OA: Solution in well (even if no Swl data we need to pump the conc)
		if (cellsHwel.size() > 0) 
		{
			Info << "     -- Head :: well boundary condition" << endl;
			for (i = 0; i < cellsHwel.size(); i++)
			{
				i1 = cellsHwel[i];
				v = mesh.V()[i1];
				unity[i1] = 0;
				if (qHwel[i] < 0)
				{
					CiSp[i1] = -qHwel[i] * (1 + 1e-6) / v; 
				} // in Sp term as discharge prop to local concentration
				else
				{
					CiSu[i1] = -qHwel[i] / v * solu_conc[condSwel[i] * ph_ncomp + ic];
				}
			}
		}

		//-- Head:: GHB boundary condition
		if (cellsHghb.size() > 0)
		{
			Info << "     -- Head :: GHB boundary condition" << endl;
			for (i = 0; i < cellsHghb.size(); i++)
			{
				i1 = cellsHghb[i];
				v = mesh.V()[i1];
				if (qHghb[i] < 0)
				{
					CiSp[i1] = -qHghb[i] / v;
				} // in Sp term as discharge prop to local concentration
				else
				{
					CiSu[i1] = -qHghb[i] / v * solu_conc[condSghb[i] * ph_ncomp + ic];
				}
			}
		} 
		
		//-- Head:: fixed head boundary condition
		if (cellsShfix.size() > 0)
		{
			Info << "     -- Head :: fixed head boundary condition" << endl;
			for (i = 0; i < cellsShfix.size(); i++)
			{
				i1 = cellsShfix[i];
				v = mesh.V()[i1]; //
				if (qHfix[i] < 0)
				{
					CiSp[i1] = -qHfix[i] / v;
				} // outflow discharge prop to local concentration
				else
				{
					CiSu[i1] = -qHfix[i] / v * solu_conc[condShfix[i] * ph_ncomp + ic];
				}
			}

		}

		//-- Head:: drain boundary condition
		//-- OA: remove conc at drn places
		if (cellsHdrn.size() > 0)
		{
			Info << "     -- Head :: drain boundary condition" << endl;
			for (i = 0; i < cellsHdrn.size(); i++)
			{
				i1 = cellsHdrn[i];
				v = mesh.V()[i1];
				CiSp[i1] = -qHdrn[i] / v;
			}
		} 
		
		//-- Head:: recharge boundary condition
		//-- OA: for recharge however there shall be a solution at all surface cells so Srch and Hrch correspond
		if (cellsHrch.size() > 0)
		{
			Info << "     -- Head :: recharge boundary condition" << endl;
			for (i = 0; i < cellsHrch.size(); i++)
			{
				i1 = clayRch[i];
				v = mesh.V()[i1];
				// recharge does not occur at drains, ghb, fix; and if rch<0: evapo, no solutes
				if ((CiSu[i1] == 0) && (qHrch[i] > 0))
				{
					CiSu[i1] += -qHrch[i] / v * solu_conc[condSrch[i] * ph_ncomp + ic];
				}
			}
		}

		//-- Do a specific treatment for dry cells
		//-- for dry cells, set conc as background
		if (flowType == 2) // unconfined flow
		{
			for (j = 0; j < ncell; j++)
			{
				if (sw[j] < sw_min[j] * 1.1)
				{
					CiSu[j] = -solu_conc[ic];
					CiSp[j] = 1;
				}
			}
		}

	}

	else // If flowType is 4, i.e. dual phase. Conditions are for water
	{
		//-- Pressure:: fixed boundary condition
		if (cellsPfix.size() > 0)
		{
			Info << "     -- Pressure :: fixed boundary condition" << endl;
			for (i = 0; i < cellsPfix.size(); i++)
			{
				i1 = cellsPfix[i];
				v = mesh.V()[i1];
				if (qPfix[i] <= 0)
				{
					CiSp[i1] = -qPfix[i] / v;
				} // outflow discharge prop to local concentration
				else
				{
					CiSu[i1] = -qPfix[i] / v * solu_conc[condSpfix[i] * ph_ncomp + ic];
				}
			}
		}

		//-- Pressure:: well boundary condition
		//- OA: solution in well (even if no Swl data we ned to pump the conc)
		if (cellsPwel.size() > 0) 
		{
			Info << "     -- Pressure :: well boundary condition" << endl;
			for (i = 0; i < cellsPwel.size(); i++)
			{
				i1 = cellsPwel[i];
				v = mesh.V()[i1];
				unity[i1] = 0;
				Info << "qPwel " << qPwel[i] << "\n";
				if (qPwel[i] <= 0)
				{
					CiSp[i1] = -qPwel[i] * (1 + 1e-6) / v;
				} // in Sp term as discharge prop to local concentration
				else
				{
					CiSu[i1] = -qFwel[i] / v * solu_conc[condSwel[i] * ph_ncomp + ic];
				}
			}
		}
	}

	//---------------------------------//
	//-- Solve the diffusion-dispersion equation for species ic
	//---------------------------------//
	Info << "   - Solve the diffusion-dispersion equation for species "<<ic<<" ... " << endl;

	//-- Concentration equation
	fvScalarMatrix CwiEqn
	(
		Rf * thetaw * fvm::ddt(Cwi) 
		+ Rf * Cwi * fvc::ddt(thetaw) 
		+ fvc::div(phiw, Cwi, "div(phiw,Cwi)") 
		- fvm::laplacian(Deffw, Cwi, "laplacian(Deffw,Cwi)")
		+ fvm::SuSp(CiSp, Cwi)
		== -CiSu
	);

	for (i = 0; i < cellsSfix.size(); i++)
	{
		CwiEqn.setReference(cellsSfix[i], solu_conc[condSfix[i] * ph_ncomp + ic]);
	}

	//-- Do a specific treatment for inactive and dry cells
	if (inactiveCells.size() + dryCells.size() > 0)
	{
		labelHashSet inactiveCellsSet(inactiveCells);
		forAll(mesh.owner(), iface)
		{
			if ((inactiveCellsSet.found(CwiEqn.lduAddr().lowerAddr()[iface])) || (inactiveCellsSet.found(CwiEqn.lduAddr().upperAddr()[iface])))
			{
				CwiEqn.upper()[iface] = 0.0;
				CwiEqn.lower()[iface] = 0.0;
			}
		}

		labelList drfList = dryFaces.shrink(); // neighb is low
		for (label faceI : drfList)
		{
			CwiEqn.lower()[faceI] = 0.0;
			CwiEqn.upper()[faceI] = 0.0;
		}
		for (label cellI : inactiveCells)
		{
			CwiEqn.diag()[cellI] = 1.0;
			CwiEqn.source()[cellI] = 0;
		}
		labelList drcList = dryCells.shrink();
		for (label cellI : drcList)
		{
			CwiEqn.diag()[cellI] = 1.0;
			CwiEqn.source()[cellI] = 0;
		}
	}

	//-- Solve the equation
	CwiEqn.solve();

	//---------------------------------//
	//-- Compute dC: deltaC
	//---------------------------------//
	Info << "   - Remove negative C and compensate mass for species "<<ic<<" ... " << endl;

	// calculating dC, and time step (only for classical species)
	if (ic > 3)
	{
		dC = 0;
		mxC = 0;
		dff = 0;
		jmax = 0;
		for (j = 0; j < nxyz; j++)
		{
			mxC = max(mxC, Cw[ic]()[ractive[j]]);
		}
		for (j = 0; j < nxyz; j++)
			if (bcCwi[j] == 0)
			{
				dff = mag(Cw[ic]()[ractive[j]] - Cw[ic]().prevIter()[ractive[j]]);
				if (dff > dC)
				{
					dC = dff;
					jmax = j;
				}
			}
		dC = dC / (mxC + SMALL);

		//-- LM:: do we need all the following output when not in debugging mode
		Info<<"     -- Overview of maximum concentration difference, max concentration and mass balance "<<endl;
		Info<<"       ---- Cell id "<<jmax<<" with |C(t)-C(t-1)|="<<dC<<endl;
		Info<<"       ---- Maximum cell concentration found at t is "<<mxC<<" kg/m3"<<endl;
		Info<<"       ---- New dt computed from concentration for this species "<<dtForC<<endl;
	
		dtForC = min(dtForC, dCmax / (max(dC, 0) + SMALL) * runTime.deltaTValue());
	}
}

//---------------------------------//
//-- Number of gas compounds
//---------------------------------//
if (ph_gcomp == 0)
{
	reactIndicator = 0.;
}

//---------------------------------//
//-- Compute the time-step
//---------------------------------//
//-- LM:: do we need all the following output when not in debugging mode
Info<< "     -- Current dt "<<runTime.deltaTValue()<< endl; 
Info<< "     -- New dt computed from flow "<<newDeltaT<< endl; 
Info<< "     -- New dt computed from concentrations "<<dtForC << endl; 
Info<< "     -- Minimum dt "<<minDeltaT << endl; 
Info<< "     -- New dt (min of flow and conc) "<<max(minDeltaT, min(dtForC, newDeltaT)) << endl;

newDeltaT = max(minDeltaT, min(dtForC, newDeltaT));

Info << "--->> End solving diffusion-dispersion equation (reactive transport - concentration Cwi). " << endl;
Info << "//******************************************************************//\n\n"
	 << endl;