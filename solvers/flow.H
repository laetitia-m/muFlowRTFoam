/**********************************************/
//-- Solve the saturation-pressure equation (Sw-p) for two-phase flow
//-- Solve the head equation (hp) for unsaturated
/**********************************************/

//---------------------------------//
//-- Variables and initialisation
//---------------------------------//
scalar CFLUse        = -1;
scalar maxDeltaTFact = -1;
scalar maxDeltaT = runTime.controlDict().lookupOrDefault("maxDeltaT", GREAT);
scalar dswmax    = runTime.controlDict().lookupOrDefault("dswmax", 0.02);

//---------------------------------//
//-- Update Uw BCs and water fluxes phiw
//---------------------------------//
Uw.correctBoundaryConditions();

forAll(mesh.boundary(), patchi)
{
	if (isA<fixedValueFvPatchField<vector>>(Uw.boundaryField()[patchi]))
	{
		phiw.boundaryFieldRef()[patchi] = Uw.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
	}
}

//---------------------------------//
//-- PICARD for unsaturated flow
//---------------------------------//
if (flowType == 3)
{
	Info << "--->> Enter PICARD for solving head equation (hp) for unsaturated..." << endl;

	//-- If the time step is controlled by Picard
	// LM: do we actually go through the other loop?
	if (timeStepControl == "Picard")
	{
		Info << "in picard " << endl;
		//-- Initialisation
		iterPicard = 0;
		resPicard  = 1000.;
		residu0    = 1000;

		// OA: This loop is imitated from Orgogozo Richards
		deltaTFact = 1.;
		sw1 = sw;

		// While the residual is higher than than tolerance set up in system/fvSolutions
		while (resPicard > tolPicard)
		{
			iterPicard++;
			sw = sw_tmp;

			//-- Update fields and solve head equation (hp)
			#include "flow/updateThetaFields.H"
			#include "flow/hpEqn.H"
			
			//-- Get the residual
			resPicard = gMax((mag(sw - sw_tmp) * 200)->internalField());

			//-- New head is relaxed
			if (residu > residu0)
			{
				hp = hp * 0.75 + hp.prevIter() * 0.25;
			}

			//-- Store the previous water saturation
			sw_tmp = sw;
			residu0 = residu;

			//-- If the max number of iterations is reached: stop.
			if (iterPicard == maxIterPicard)
			{
				runTime.setTime(presentTime, tstep);
				break;
			}
		}

		Info << "Picard nb iterations : " << iterPicard << endl;
		if (iterPicard == maxIterPicard)
		{
			deltaTFact = 0.2;
			sw_tmp = sw1;
			sw = sw1;
			rewind = 1;
		} // hp_tmp=hp1;hp=hp1;rewind=1;
		else if (iterPicard <= nIterStability)
		{
			deltaTFact = 1.2;
		}
		newDeltaT = min(deltaTFact * newDeltaT, maxDeltaT); // runTime.deltaTValue()
		#include "flow/updateThetaFields.H"
	}
	else
	{
		#include "flow/hpEqn.H"
		#include "flow/updateThetaFields.H"
		volScalarField::Internal dsw_tmp = mag(sw - sw.oldTime());
		dsw = gMax(dsw_tmp) / average(sw).value(); // cannot do directly linne aboe and that one
		scalar ratio = dswmax / (dsw + SMALL);
		Info << "dsw " << dsw << " dswmax " << dswmax << " ratio " << ratio << endl;
		deltaTFact = max(min(ratio, 1.2), 0.1);
		newDeltaT = min(deltaTFact * runTime.deltaTValue(), maxDeltaT);
	}
}

//---------------------------------//
//-- PICARD for two-phase flow
//---------------------------------//
else if (flowType == 4)
{
	Info << "--->> Enter PICARD for solving saturation-pressure equation (S-p) for 2-phase flow..." << endl;

	if (timeStepControl == "Picard")
	{
		Info << "in picard " << endl;
		iterPicard = 0;
		resPicard = 1000.;
		// this loop is imitated from Orgogozo Richards
		deltaTFact = 1.;
		p1 = p * 1;
		while (resPicard > tolPicard)
		{
			iterPicard++;
			p_tmp = p;
			#include "flow/pSEqn.H"
			resPicard = gMax((mag(p - p_tmp) / atmPa * 1e3)->internalField());
			Info << "res Pic " << resPicard << endl;
			if (iterPicard == maxIterPicard)
			{
				break;
			}
		} // end picard iter
		Info << "Picard nb iterations : " << iterPicard << endl;
		if (iterPicard == maxIterPicard)
		{
			deltaTFact = 0.2;
			p = p1;
			rewind = 1;
		}
		else if (iterPicard <= nIterStability)
		{
			deltaTFact = 1.2;
		}
		newDeltaT = min(deltaTFact * runTime.deltaTValue(), maxDeltaT);
	}

	else
	{
		#include "flow/pSEqn.H"
		volScalarField::Internal dsw_tmp = mag(sw - sw.oldTime());
		for (j = 0; j < nxyz; j++)
		{
			if (bcPS[j] == 0)
			{
				dsw_tmp[j] = 0;
			}
		}
		dsw = gMax(dsw_tmp) / average(sw).value(); // cannot do directly linne aboe and that one
		scalar ratio = dswmax / (dsw + SMALL);
		deltaTFact = max(min(ratio, 1.2), 0.1);
		newDeltaT = min(deltaTFact * runTime.deltaTValue(), newDeltaT);
	}
	Info << " dswmax " << dswmax << " ndetlaT " << newDeltaT << endl;
}

//---------------------------------//
//-- PICARD for no flow
//---------------------------------//
else if (flowType != 0) // flowType=0 no flow
{
	Info << "--->> Enter PICARD for solving head equation (h) for no flow..." << endl;


	Info << "in simple flow " << endl;
	#include "flow/hEqn.H"
	newDeltaT = min(maxDeltaT, 1.4 * newDeltaT);
}

Info << " deltaTFact " << deltaTFact << " maxDeltaT " << maxDeltaT << " newdeltat " << newDeltaT << endl;

Uw = fvc::reconstruct(phiw); // This is a true darcy velocity (flux divided by surface even in unconfined)
