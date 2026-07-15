/**********************************************/
//-- Solve for head equation - variable hp
//-- For unsaturated
/**********************************************/

Info << "\n\n//*******************************************************************//" << endl;
Info << "//* * * * * SOLVE HEAD EQUATION (hp) for unsaturated * * * * *//" << endl;
Info << "//*******************************************************************//" << endl;
Info << "--->> Start solving head equation (hp) for unsaturated..." << endl;

//---------------------------------//
//-- Variables and initialisation
//---------------------------------//

//-- Time
time = mesh.time().value(); 

//-- Source terms
hSu = hSu * 0;
hSp = hSp * 0;

//-- Store previous solution
hp.storePrevIter();
sw.storePrevIter();

//---------------------------------//
//-- Boundary conditions
//---------------------------------//
Info << "   - Set up the boundary conditions... " << endl;

//-- Head well boundary condition
if (cellsHwel.size() > 0)
{
	Info << "     -- Head:: well boundary condition" << endl;
	if (time >= tHwel)
	{
		output = getCbuffer("hwel", itHwel, cellsHwel.size());
		tHwel = output.t;
		cellsHwelData = output.d;
		itHwel += 1;
		tnext = min(tnext, tHwel);
	}
	for (size_t i = 0; i < cellsHwel.size(); i++)
	{
		i1 = cellsHwel[i];
		qHwel[i] = cellsHwelData[i * 4 + 2];
		hSu[i1] = qHwel[i] / mesh.V()[i1];
	}
}

//-- Drain boundary condition
if (cellsHdrn.size() > 0)
{
	Info << "     -- Head:: drain boundary condition" << endl;
	if (time >= tHdrn)
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
		v = mesh.V()[i1];
		hSu[i1] = cellsHdrnData[i * 4 + 2] / v;
		hSp[i1] = cellsHdrnData[i * 4 + 3] / v;
		if (hp()[i1] < cellsHdrnData[i * 4 + 2] / cellsHdrnData[i * 4 + 3])
		{
			hSu[i1] /= 1e9;
			hSp[i1] /= 1e9;
		} // almost no flow if head below fixed one
	}
}

//-- Recharge boundary condition
if (cellsHrch.size() > 0)
{
	Info << "     -- Head:: recharge boundary condition" << endl;
	if (time >= tHrch)
	{
		output = getCbuffer("hrch", itHrch, cellsHrch.size());
		tHrch = output.t;
		cellsHrchData = output.d;
		itHrch += 1;
		tnext = min(tnext, tHrch);
	}

	//-- Local variables
	float frech = 1; // frech is used to stop recharge or evapo when hp is too high or too low
	float hplim = -1.0 / lunits;
	for (size_t i = 0; i < cellsHrch.size(); i++)
	{
		i1 = cellsHrch[i];
		clayRch[i] = i1;
		qHrch[i] = cellsHrchData[i * 4 + 2]; // The injection is at the exact place of recharge (not like in unconfined)
		
		//-- Apply recharge only if there is nothing else in the cell
		if (hSu[i1] == 0)
		{
			if (qHrch[i] > 0)
			{
				frech = 1;
			}
			else
			{
				if (hp()[i1] < hplim) // hp<-3 m stop all evt
				{
					frech = 0;
				}
				else if (hp()[i1] < 0) // correction to limit head, if cm, lunits= 0.01 -> 40 cm head
				{
					frech = 0.5 * (1 + std::erf((hp()[i1] + 0.4 / lunits) / 20));
				} 
				else
				{
					frech = 1;
				}
			}
			hSu[i1] = frech * qHrch[i] / mesh.V()[i1];
		}
	}
}

//---------------------------------//
//-- Solve the head equation hp
//---------------------------------//
Info << "   - Solve the head equation (hp)... " << endl;

//-- Head equation
fvScalarMatrix hpEqn
(
	stor *fvm::ddt(hp) 
	- fvm::laplacian(Mf, hp) 
	+ fvc::div(phiGr) 
	+ fvm::SuSp(hSp, hp) 
	== hSu
);

//-- Get the fixed head conditions
// LM:: not clear why this is needed
if (cellsHfix.size() > 0)
{
	if (time >= tHfix)
	{
		output = getCbuffer("hfix", itHfix, cellsHfix.size());
		tHfix = output.t;
		cellsHfixData = output.d;
		itHfix += 1;
		tnext = min(tnext, tHrch);
	}
	for (size_t i = 0; i < cellsHfix.size(); i++)
	{
		hpEqn.setReference(cellsHfix[i], cellsHfixData[i * 4 + 2]);
	}
}

//-- Solve and Residual
residu = hpEqn.solve().initialResidual(); // Solve

Info << "     -- Resiual:: " << residu
		<< ", Execution Time:: " << runTime.elapsedCpuTime() << " s"
		<< endl;

//---------------------------------//
//-- Compute water fluxes and velocity
//---------------------------------//
Info << "   - Water fluxes and velocity computation... \n" << endl;
phiw = phiGr - (Mf * fvc::snGrad(hp)) * mesh.magSf();
Uw = fvc::reconstruct(phiw);

//-- Recompute fluxes at boundaries for fixed and drain BCs
volScalarField q1 = fvc::div(phiw);

for (size_t i = 0; i < cellsHfix.size(); i++)
{
	i1 = cellsHfix[i];
	qHfix[i] = q1[i1] * mesh.V()[i1];
} 

for (size_t i = 0; i < cellsHdrn.size(); i++)
{
	i1 = cellsHdrn[i];
	qHdrn[i] = (hSu[i1] - hp()[i1] * hSp[i1]) * mesh.V()[i1];
}

Info << "--->> End solving head equation (h) for unsaturated. " << endl;
Info << "//******************************************************************//\n\n" << endl;
