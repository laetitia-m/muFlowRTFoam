/**********************************************/
//-- Solve for head equation - variable h
//-- For confined and unconfined
/**********************************************/

Info << "\n\n//******************************************************************//" << endl;
Info << "//* * * SOLVE HEAD EQUATION (h) for confined and/or unconfined * * *//" << endl;
Info << "//******************************************************************//" << endl;
Info << "--->> Start solving head equation (h) for confined and/or unconfined..." << endl;

/*
Here we calculate the head equation for different conditions (confined, unconfined) the unsaturated is calculated in hpEqn,
using the boundary conditions
the fluxes at boundaries are calculated in the budget file
*/

//---------------------------------//
//-- Variables and initialisation
//---------------------------------//
double x; // relative saturated thickness
int ilay, nc;
double e = 1e-5;

//-- Time
time = mesh.time().value();

//-- Source terms
hSu = hSu * 0;
hSp = hSp * 0;

//-- Dry cells
dryCells.clear();
dryFaces.clear();

//-- Store previous solution
h.storePrevIter();

//---------------------------------//
//-- Compute the relative thickness 
//-- water saturant and mobility
//-- for unconfined flow
//---------------------------------//
Mf = M0f * 1;

if (flowType == 2) // unconfined
{
	//-- Find Sw
	for (ilay = 0; ilay < nlay; ilay++)
	{
		for (int i = 0; i < ncell_lay; i++)
		{
			nc = ilay * ncell_lay + i;
			x = (h[nc] - zbot[nc]) / thk[nc]; // relative thickness
			if (x < 0)
			{
				dryCells.append(nc);
			}
			if (x < e && x > 0)
			{
				sw[nc] = pow(x, 2) / (4 * e) + x / 2 + e / 4; // OA: same eq as in mdf USg (transport manual fg 21). LM: ??
			} 
			else if (x > 1 - e && x < 1)
			{
				sw[nc] = 1 - pow((1 - x), 2) / (4 * e) - (1 - x) / 2 - e / 4;
			}
			else
			{
				sw[nc] = std::max(std::min(x, sw_max[nc]), sw_min[nc]);
			}
		}
	}

	//-- Compute mobility Mf
	forAll(mesh.owner(), iface) // calculate Mf
	{
		oCell = mesh.owner()[iface];
		odCface = mag(oCf[iface] - oC[oCell]);
		nCell = mesh.neighbour()[iface];
		ndCface = mag(oCf[iface] - oC[nCell]);

		//-- If owner is dry or neighbour is dry
		if ((sw[oCell] < sw_min[oCell] * 1.1) || (sw[nCell] < sw_min[nCell] * 1.1))
		{
			dryFaces.append(iface);
		}															  // own or nghb is dry
		if (abs(mesh.Sf()[iface].z()) < mag(mesh.Sf()[iface]) * 1e-5) // vertical faces
		{
			scalar swm = (sw[oCell] * sw[nCell]) / (odCface * sw[oCell] + ndCface * sw[nCell]) * (ndCface + odCface);
			Mf[iface] = M0f[iface] * swm;
		}
		else // other faces if sw=sw_min Mf=Small
		{
			if ((sw[oCell] < sw_min[oCell] * 1.1) || (sw[nCell] < sw_min[nCell] * 1.1))
			{
				Mf[iface] = M0f[iface] * sw_min[oCell];
			}
		}
	}
}

//---------------------------------//
//-- Compute the storage coefficient
//---------------------------------//
forAll(stor(), celli)
{
	if (sw[celli] < sw_max[nc])
	{
		stor[celli] = storU[celli] / thk[celli];
	}
	else
	{
		stor[celli] = storC[celli] / thk[celli];
	}
};
stor.dimensions().reset(dimless / dimLength);

//---------------------------------//
//-- Boundary conditions
//---------------------------------//
Info << "   - Set up the boundary conditions... " << endl;

//-- Head well boundary condition
if (cellsHwel.size() > 0)
{
	Info << "     -- Head:: well boundary condition" << endl;
	if (time > tHwel)
	{
		output = getCbuffer("hwel", itHwel, cellsHwel.size());
		tHwel = output.t;
		cellsHwelData = output.d;
		itHwel += 1;
		tnext = min(tnext, tHwel);
	}
	for (size_t i = 0; i < cellsHwel.size(); i++)
	{
		v = mesh.V()[cellsHwel[i]];
		qHwel[i] = cellsHwelData[i * 4 + 2];
		hSu[cellsHwel[i]] = qHwel[i] / v;
	}
}

//-- Head GHB boundary condition
if (cellsHghb.size() > 0)
{
	Info << "     -- Head:: GHB boundary condition" << endl;
	if (time > tHghb)
	{
		output = getCbuffer("hghb", itHghb, cellsHghb.size());
		tHghb = output.t;
		cellsHghbData = output.d;
		itHghb += 1;
		tnext = min(tnext, tHghb);
	}
	for (size_t i = 0; i < cellsHghb.size(); i++)
	{
		i1 = cellsHghb[i];
		v = mesh.V()[i1];
		hSu[i1] = cellsHghbData[i * 4 + 2] / v;
		hSp[i1] = cellsHghbData[i * 4 + 3] / v;
	}
}

//-- Head river boundary condition
if (cellsHriv.size() > 0)
{
	Info << "     -- Head:: river boundary condition" << endl;
	if (time > tHriv)
	{
		output = getCbuffer("hriv", itHriv, cellsHriv.size());
		tHriv = output.t;
		cellsHrivData = output.d;
		itHriv += 1;
		tnext = min(tnext, tHriv);
	}
	for (size_t i = 0; i < cellsHriv.size(); i++)
	{
		i1 = cellsHriv[i];
		v = mesh.V()[i1];
		hSu[i1] = cellsHrivData[i * 4 + 2] / v;
		hSp[i1] = cellsHrivData[i * 4 + 3] / v;
	}
}

//-- Head drain boundary condition
if (cellsHdrn.size() > 0)
{
	Info << "     -- Head:: drain boundary condition" << endl;
	if (time > tHdrn)
	{
		output = getCbuffer("hdrn", itHdrn, cellsHdrn.size());
		tHdrn = output.t;
		cellsHdrnData = output.d;
		itHdrn += 1;
		tnext = min(tnext, tHdrn);
	}
	for (size_t i = 0; i < cellsHdrn.size(); i++)
	{
		i1 = cellsHdrn[i];
		v = mesh.V()[i1]; // missing condition!!!
		hSu[i1] = cellsHdrnData[i * 4 + 2] / v;
		hSp[i1] = cellsHdrnData[i * 4 + 3] / v;
		if (h()[i1] < cellsHdrnData[i * 4 + 2] / cellsHdrnData[i * 4 + 3])
		{
			hSu[i1] /= 1e6;
			hSp[i1] /= 1e6;
		} // no flow if head below fixed one
	}
}

//-- Head recharge boundary condition
if (cellsHrch.size() > 0)
{
	Info << "     -- Head:: recharge boundary condition" << endl;
	if (time > tHrch)
	{
		output = getCbuffer("hrch", itHrch, cellsHrch.size());
		tHrch = output.t;
		cellsHrchData = output.d;
		itHrch += 1;
		tnext = min(tnext, tHrch);
	}
	int ibl, Rblok, nRblok;
	if (ncell_row > 1)
	{
		Rblok = ncell_row;
		nRblok = ncell_col;
	} // case xsection
	else
	{
		Rblok = ncell_lay;
		nRblok = nlay;
	} // 3d case
	for (size_t i = 0; i < cellsHrch.size(); i++)
	{
		i1 = cellsHrch[i];
		for (ibl = 0; ibl < nRblok; ibl++)
		{
			nc = i1 - ibl * Rblok;
			if (sw[nc] >= sw_max[nc])
			{
				break;
			};
		}
		
		//-- Search for the first confined layer from the top
		int ibl1 = max(0, ibl - 1); // ibl1 from top
		nc = i1 - ibl1 * Rblok;
		clayRch[i] = nc;
		qHrch[i] = cellsHrchData[i * 4 + 2];
		if (hSu[nc] == 0)
		{
			hSu[nc] = qHrch[i] / mesh.V()[nc];
		} // injection of recharge in the lowest unconfined layer only if no fix cond here
	}
}

//-- Head fixed head boundary condition
if (cellsHfix.size() > 0)
{
	Info << "     -- Head:: fixed head boundary condition" << endl;
	if (time > tHfix)
	{
		output = getCbuffer("hfix", itHfix, cellsHfix.size());
		tHfix = output.t;
		cellsHfixData = output.d;
		itHfix += 1;
		tnext = min(tnext, tHfix);
	}
	for (size_t i = 0; i < cellsHfix.size(); i++)
	{
		i1 = cellsHfix[i];
		v = mesh.V()[i1];
		hSu[i1] = cellsHfixData[i * 4 + 2] * 1e3;
		hSp[i1] = 1e3; // high C conductance
	}
}

//---------------------------------//
//-- For coupling, compute the multiplying factors
//---------------------------------//
//-- Get the multiplying factors from plugin_H
if (coupling == 1)
{
	volScalarField factM = plugH.corrM(mesh, freak);
	surfaceScalarField factMf = fvc::interpolate(factM, "factM");
	Mf *= factMf;

	//-- LM: the followinv seems to be only terminal output? Do we really need the following?
	for (int i = 0; i < 6; i++)
	{
		Info << factM[i] << " ";
	}
	if (flagW == 1)
	{
		std::cout << "in write cpl fact ??? " << "\n";
		std::ofstream outFact(cur_dir / runTime.timeName() / "factM");
		outFact.unsetf(std::ios::scientific);
		outFact.precision(6);
		for (int j = 0; j < ncell; j++)
		{
			outFact << factM[j] << "\n";
		}
		outFact << std::endl;
		outFact.close();
	}
}

//---------------------------------//
//-- Solve the head equation h
//---------------------------------//
Info << "   - Solve the head equation (h)... " << endl;

//-- Head equation
fvScalarMatrix hEqn
(
	stor *fvm::ddt(h) 
	- fvm::laplacian(Mf, h) 
	+ fvm::SuSp(hSp, h) 
	== hSu
);

//-- Solve and Residual
residu = hEqn.solve().initialResidual(); // Solve

Info << "     -- Resiual:: " << residu
		<< ", Execution Time:: " << runTime.elapsedCpuTime() << " s"
		<< endl;

//---------------------------------//
//-- Compute water fluxes, velocity and BCs
//---------------------------------//
Info << "   - Water fluxes and velocity computation... \n" << endl;
phiw == (-Mf * fvc::snGrad(h)) * mesh.magSf(); // water fluxes
Uw = fvc::reconstruct(phiw); // water velocity

//-- Recompute fluxes at boundaries for fixed, drain and GHB BCs
volScalarField q1 = fvc::div(phiw);

for (size_t i = 0; i < cellsHfix.size(); i++)
{
	i1 = cellsHfix[i];
	qHfix[i] = q1[i1] * mesh.V()[i1];
} 

for (size_t i = 0; i < cellsHdrn.size(); i++)
{
	i1 = cellsHdrn[i];
	qHdrn[i] = (hSu[i1] - h()[i1] * hSp[i1]) * mesh.V()[i1];
}
for (size_t i = 0; i < cellsHghb.size(); i++)
{
	i1 = cellsHghb[i];
	qHghb[i] = (hSu[i1] - h()[i1] * hSp[i1]) * mesh.V()[i1];
}

Info << "--->> End solving head equation (h) for confined and/or unconfined. " << endl;
Info << "//******************************************************************//\n\n" << endl;