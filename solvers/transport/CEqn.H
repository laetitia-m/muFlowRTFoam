/**********************************************/
//-- Solve the diffusion-dispersion equation
//-- For concentration C when no reaction
/**********************************************/

Info << "\n\n//******************************************************************//" << endl;
Info << "//* * * * * * SOLVE DIFFUSION-DISPERSION EQUATION C * * * * * * * * //" << endl;
Info << "//******************************************************************//" << endl;
Info << "--->> Start solving diffusion-dispersion equation (no reaction - concentration C) " << endl;

//---------------------------------//
//-- Variables and initialisation
//---------------------------------//
size_t i;
dimensionedScalar USMALL("",dimLength/dimTime,1e-36); // To avoid dividing by O in effective water diffusivity
double dff=0;
int imax=0;

//-- Initialisation
a=a*0;
dC=0;
mxC=0;
sC=0;
sC1=0;

//- Time
float time = mesh.time().value();

//-- Source terms
CSu=CSu*0;
CSp=CSp*0;

//-- Water content
thetaw=eps*sw;

//-- Store previous solution
C.storePrevIter();


//---------------------------------//
//-- Compute the effective water diffusivity
//---------------------------------//
Deffw = (alphaT * mag(Uw) + Dw0* pow(thetaw,3.333) /pow(eps,2.))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL);
Info << "   - Diffusivity: " << Dw0.value() << endl;
Info << "   - Effective diffusivity matrix: " << Deffw[15] << endl;

//---------------------------------//
//-- Boundary conditions
//---------------------------------//
Info << "   - Set up the boundary conditions... " << endl;

//-- Head well boundary condition
if (cellsHwel.size()>0)
{
	Info << "     -- Head:: well boundary condition" << endl;
	if (time>tCwel)
	{
		output = getCbuffer("cwel",itCwel,cellsCwel.size()); 
		tCwel  = output.t;
		cellsCwelData=output.d;
		itCwel +=1;
		tnext   = min(tnext,tCwel);
	}
	for (i=0; i<cellsHwel.size(); i++) 
	{ 
		i1=cellsCwel[i];
		v=mesh.V()[i1];
		if (qHwel[i]<0) 
		{
			CSp[i1] = -qHwel[i]/v;
		} // In Sp term as discharge prop to local concentration, a positive Sp will make solution implicit
		else 
		{
			CSu[i1] = -qHwel[i]/v*cellsCwelData[i*4+2];
		} // The injected mass is the water flow x the conc
	}
}

//-- Head GHB boundary condition
//-- OA: Hghb because Cghb may not be set here
if (cellsHghb.size()>0)
{
	Info << "     -- Head:: GHB boundary condition" << endl;
	if (time>tCghb) 
	{
		output = getCbuffer("cghb",itCghb,cellsCghb.size()); 
		tCghb  = output.t;
		cellsCghbData=output.d;
		itCghb +=1;
		tnext = min(tnext,tCghb);
	}
	for (i=0; i<cellsHghb.size(); i++) 
	{ 
		i1=cellsHghb[i];
		v=mesh.V()[i1];
		if (qHghb[i]<0) 
		{
			CSp[i1] = -qHghb[i]/v;
		} // outflow discharge prop to local concentration
		else 
		{
			CSu[i1] = -qHghb[i]/v*cellsCghbData[i*4+2];
		}
	}
}

//-- Head fixed head boundary condition
//-- OA: store the chfix conditions data
if (cellsHfix.size()>0)
{
	Info << "     -- Head:: fixed head boundary condition" << endl;
	if (time>tChfix) 
	{
		output = getCbuffer("chfix",itChfix,cellsChfix.size()); 
		tChfix = output.t;
		cellsChfixData=output.d;
		itChfix +=1;
		tnext = min(tnext,tChfix);
	}
	for (i=0; i<cellsHfix.size(); i++) 
	{ 
		i1=cellsHfix[i];
		v=mesh.V()[i1];
		if (qHfix[i]<0) 
		{
			CSp[i1] = -qHfix[i]/v;
		} // outflow discharge prop to local concentration
		else 
		{
			CSu[i1] = -qHfix[i]/v*cellsChfixData[i*4+2];
		}
	}
}

//-- Head drain boundary condition
//-- OA: for drain cells (only out), same as above, except that drains are constant in time
if (cellsHdrn.size()>0)
{
	Info << "     -- Head:: drain boundary condition" << endl;
	for (i=0; i<cellsHdrn.size(); i++) 
	{ 
		i1=cellsHdrn[i];
		v=mesh.V()[i1];
		if (qHdrn[i]<0) 
		{
			CSp[i1] = -qHdrn[i]/v;
		} // outflow qHdrn<0 discharge prop to local concentration
	}
}

//-- Concentration recharge boundary condition
//-- OA: recharge is different, required only if C>0
if (cellsCrch.size()>0) 
{
	Info << "     -- Concentration:: recharge boundary condition" << endl;
	if (time>tCrch) 
	{
		output = getCbuffer("crch",itCrch,cellsCrch.size()); 
		tCrch = output.t;
		cellsCrchData=output.d;
		itCrch +=1;
		tnext = min(tnext,tCrch);
	}
	for (i=0; i<cellsHrch.size(); i++) 
	{ 
		i1=clayRch[i];
		v=mesh.V()[i1];
		if ((CSu[i1]==0)&&(qHrch[i]>0)) //recharge occurs only where no ghb, riv or drains + rch<0=evapo and no solutes evporaiton
		{ 
			CSu[i1] = -qHrch[i]/v*cellsCrchData[i*4+2];
		}
	} 
}

//-- Fixed concentration boundary condition
//-- OA: fixed cells done with ~ghb with very high conductance
if (cellsCfix.size()>0)
{
	Info << "     -- Concentration:: Fixed concentration boundary condition" << endl;
	if (time>tCfix) 
	{
		output = getCbuffer("cfix",itCfix,cellsCfix.size()); 
		tCfix = output.t;
		cellsCfixData=output.d;
		itCfix +=1;
		tnext = min(tnext,tCfix);
	}
	for (i=0; i<cellsCfix.size(); i++) 
	{ 
		i1=cellsCfix[i];
		v=mesh.V()[i1];
		CSu[i1] = -cellsCfixData[i*4+2]*1e5;
		CSp[i1] = 1e5;// high C conductance
	}
}
	
//---------------------------------//
//-- Solve the diffusion-dispersion equation
//---------------------------------//
Info << "   - Solve the diffusion-dispersion equation for C... " << endl;

//-- Concentration equation
fvScalarMatrix CEqn
(
	Rf * thetaw * fvm::ddt(C) 
	+ Rf * C * fvc::ddt(thetaw) 
	+ fvc::div(phiw, C, "div(phiw,C)") 
	- fvm::laplacian(Deffw, C, "laplacian(Deffw,C)")
	+ thetaw * fvm::Sp(lbdaw,C)
	+ fvm::SuSp(CSp, C)
	== - CSu - thetaw * gamw
);

//-- Do a specific treatment for inactive and dry cells
if (inactiveCells.size()+dryCells.size()>0) 
{
	Info<<"     -- In inactive/dry cells:: nbr inactive cells "<<inactiveCells.size()<<" & dry cells "<<dryCells.size()<<endl;
	labelHashSet inactiveCellsSet(inactiveCells);
	labelList drcList = dryCells.shrink();
	
	//-- OA: ChatGPT
	forAll(dryFaces, i)
	{
		label faceI = dryFaces[i];
		label owner = mesh.owner()[faceI];
		label neigh = mesh.neighbour()[faceI];

		//-- Enlever contribution diagonale AVANT de couper
		CEqn.diag()[owner] -= CEqn.upper()[faceI];
		CEqn.diag()[neigh] -= CEqn.lower()[faceI];

		//-- Couper le lien
		CEqn.upper()[faceI] = 0.0;
		CEqn.lower()[faceI] = 0.0;
	}

	for (label cellI : drcList) 
	{
		CEqn.diag()[cellI] = 1.0;
		CEqn.source()[cellI] = 0.0;
	}
}

//-- Solve the equation
CEqn.solve();

//---------------------------------//
//-- Compute dC: deltaC
//---------------------------------//
//-- Removing negative value increase artifially the total
Info << "   - Remove negative C and compensate mass... " << endl;

//-- Removing negative C and compensating mass
for (int i=0;i<ncell;i++) 
{
	mxC=max(mxC,C()[i]);
	sC+=C()[i]*eps[i]*mesh.V()[i];
}

C=max(C,a); // removing negative

for (int i=0;i<ncell;i++)
{
	dff = mag(C()[i]-C.prevIter()[i]);
	if (dff>dC)
	{
		dC=dff;
		imax=i;
	}
	sC1+=C()[i]*eps[i]*mesh.V()[i];
}

//-- LM:: do we need all the following output when not in debugging mode
Info<<"     -- Overview of maximum concentration difference, max concentration and mass balance "<<endl;
Info<<"       ---- Cell id "<<imax<<" with |C(t)-C(t-1)|="<<dC<<endl;
Info<<"       ---- Maximum cell concentration found at t is "<<mxC<<" kg/m3"<<endl;
Info<<"       ---- Mass computed with    negative value "<<sC<<endl;
Info<<"       ---- Mass computed without negative value "<<sC1<<endl;
Info<<"       ---- Ratio of mass "<< sC/(sC1+SMALL)<<" & Delta mass "<< sC-sC1<<endl;

//- Ratio should be equal to 1 if no mass added artificially
C *= sC/(sC1+SMALL);
dC = dC/(mxC+SMALL);

//---------------------------------//
//-- Compute the time-step
//---------------------------------//
Info << "   - Compute the time step... " << endl;

dtForC = dCmax/(max(dC,0)+SMALL)*runTime.deltaTValue(); // LM: need to better understand this condition

//-- LM:: do we need all the following output when not in debugging mode
Info<< "     -- Current dt "<<runTime.deltaTValue()<< endl; 
Info<< "     -- New dt computed from flow "<<newDeltaT<< endl; 
Info<< "     -- New dt computed from concentration "<<dtForC << endl; 
Info<< "     -- New dt (min of flow and conc) "<<min(dtForC, newDeltaT) << endl;

newDeltaT = min(dtForC, newDeltaT);

Info << "--->> End solving diffusion-dispersion equation (no reaction - concentration C). " << endl;
Info << "//******************************************************************//\n\n" << endl;