/**********************************************/
//-- Solve for steady head equation - variable h
/**********************************************/

Info << "\n\n//******************************************************************//" << endl;
Info << "//* * * * * * * * SOLVE STEADY HEAD EQUATION (h) * * * * * * * * * *//" << endl;
Info << "//******************************************************************//" << endl;
Info << "--->> Start solving steady head equation (h)..." << endl;

//---------------------------------//
//-- Variables
//---------------------------------//
simpleControl simple(mesh);

double x; // relative saturated thickness
double e = 1e-5;
int ilay, nc;

residu = 1e5;
surfaceScalarField M1f = M0f;

for (int ii = 0; ii < 50; ii++) // LM: why 50?
{
	// a: relaxation factor?
	double a = std::min(double(ii / 20.), 1.); // LM: a will always be ii/20.
	Info << "   - Iteration: " << ii << endl;

	/* if they are bc in Uw :-> removed from now one
	Uw.correctBoundaryConditions();
	forAll(mesh.boundary(),patchi)
		{
		if (isA< fixedValueFvPatchField<vector> >(Uw.boundaryField()[patchi]))
			{ phiw.boundaryFieldRef()[patchi] = Uw.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi]; }
		}
	*/

	//---------------------------------//
	//-- For unconfined flow
	//---------------------------------//
	//-- Compute the effective thickness and mobility
	if (flowType == 2)
	{
		Info << "     -- For unconfined flow, compute relative saturated thickness and mobilities... " << endl;
		for (ilay = 0; ilay < nlay; ilay++) // ilay will be used from top
		{
			for (i = 0; i < ncell_lay; i++)
			{
				nc = ilay * ncell_lay + i;
				x = (h[nc] - zbot[nc]) / thk[nc]; // Relative saturated thickness
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

		forAll(mesh.owner(), iface)
		{
			oCell = mesh.owner()[iface];
			odCface = mag(oCf[iface] - oC[oCell]);
			nCell = mesh.neighbour()[iface];
			ndCface = mag(oCf[iface] - oC[nCell]);
			if (abs(mesh.Sf()[iface].z()) < mag(mesh.Sf()[iface]) * 1e-5) // Vertical faces
			{
				scalar swm = (sw[oCell] * sw[nCell]) / (odCface * sw[oCell] + ndCface * sw[nCell]) * (ndCface + odCface);
				Mf[iface] = M0f[iface] * swm;
			}
			else
			{
				Mf[iface] = M0f[iface];
			}
		}
	}

	//---------------------------------//
	//-- Compute the mobility
	//---------------------------------//
	//-- Use of relaxation factor a
	M1f = M0f * (1 - a) + Mf * a;

	//---------------------------------//
	//-- Boundary conditions
	//---------------------------------//
	//-- OA: (can be simplified only the 1st time is needed)
	//-- Initialisation
	hSu = hSu * 0;
	hSp = hSp * 0;
	size_t i; // previous cells with values must be destroyed

	Info << "     -- Set up the boundary conditions... " << endl;

	//-- Well discharge boundary condition
	if (cellsHwel.size() > 0)
	{
		Info << "       ---- Well discharge boundary condition. " << endl;
		for (i = 0; i < cellsHwel.size(); i++)
		{
			i1 = cellsHwel[i];
			v = mesh.V()[i1];
			hSu[i1] = cellsHwelData[i * 4 + 2] / v;
		}
	}

	//-- GHB (general head) condition
	if (cellsHghb.size() > 0)
	{
		Info << "       ---- GHB boundary condition " << endl;
		for (i = 0; i < cellsHghb.size(); i++)
		{
			i1 = cellsHghb[i];
			v = mesh.V()[i1];
			hSu[i1] = cellsHghbData[i * 4 + 2] / v;
			hSp[i1] = cellsHghbData[i * 4 + 3] / v;
		}
	}

	//-- River boundary condition
	if (cellsHriv.size() > 0)
	{
		Info << "       ---- River boundary condition " << endl;
		for (i = 0; i < cellsHriv.size(); i++)
		{
			i1 = cellsHriv[i];
			v = mesh.V()[i1];
			hSu[i1] = cellsHrivData[i * 4 + 2] / v;
			hSp[i1] = cellsHrivData[i * 4 + 3] / v;
		}
	}

	//-- Drain boundary condition
	if (cellsHdrn.size() > 0)
	{
		Info << "       ---- Drain boundary condition " << endl;
		for (i = 0; i < cellsHdrn.size(); i++)
		{
			i1 = cellsHdrn[i];
			v = mesh.V()[i1];
			hSu[i1] = cellsHdrnData[i * 4 + 2] / v;
			hSp[i1] = cellsHdrnData[i * 4 + 3] / v;

			//-- No flow if head below fixed one
			if (h()[i1] < cellsHdrnData[i * 4 + 2] / cellsHdrnData[i * 4 + 3])
			{
				hSu[i1] /= 1000.;
				hSp[i1] /= 1000.;
			}
		}
	}

	//-- Recharge boundary condition
	if (cellsHrch.size() > 0)
	{
		Info << "       ---- Recharge boundary condition " << endl;

		int ibl, Rblok, nRblok;
		if (ncell_row > 1) // case xsection
		{
			Rblok = ncell_row;
			nRblok = ncell_col;
		}
		else // 3d case
		{
			Rblok = ncell_lay;
			nRblok = nlay;
		}

		for (i = 0; i < cellsHrch.size(); i++)
		{
			i1 = cellsHrch[i];
			//-- Search for the first confined layer from the top
			for (ibl = 0; ibl < nRblok; ibl++)
			{
				nc = i1 - ibl * Rblok;
				if (sw[nc] >= sw_max[nc])
				{
					break;
				};
			}
			int ibl1 = max(0, ibl - 1); // ibl1 from top
			nc = i1 - ibl1 * Rblok;
			if (hSu[nc] == 0)
			{
				hSu[nc] += cellsHrchData[i * 4 + 2] / mesh.V()[nc]; // injeciton of recharge in the lowest unconfined layer
			}
		}
	}

	//---------------------------------//
	//-- Solve the steady head equation
	//---------------------------------//
	Info << "     -- Solve the steady head equation (h)... " << endl;

	//-- Steady head equation
	fvScalarMatrix hEqnStd
	(
		-fvm::laplacian(Mf, h) 
		+ fvm::Sp(hSp, h) 
		== hSu
	);

	//-- Get the fixed head conditions
	if (cellsHfix.size() > 0)
	{
		for (i = 0; i < cellsHfix.size(); i++)
		{
			hEqnStd.setReference(cellsHfix[i], cellsHfixData[i * 4 + 2]);
		}
	}

	//-- Solve and Residual
	residu = hEqnStd.solve().initialResidual(); // Solve
	Info << "     -- Resiual:: " << residu
		 << ", Relaxation factor:: " << a
		 << ", Execution Time:: " << runTime.elapsedCpuTime() << " s"
		 << endl;

	//-- Stop condition: residual<1e-6 and at least 10 iteration
	if ((ii > 10) && (residu < 1e-6))
	{
		break;
	}

	//---------------------------------//
	//-- Compute water fluxes and velocity
	//---------------------------------//
	Info << "     -- Water fluxes and velocity computation... \n" << endl;
	phiw = hEqnStd.flux();
	Uw = fvc::reconstruct(phiw);
}

Info << "---<< End of solving steady head equation (h)." << endl;
Info << "//******************************************************************//\n\n"<< endl;
