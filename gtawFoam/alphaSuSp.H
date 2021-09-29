boundedAlpha1 = min(max(alpha1, scalar(0)), scalar(1));
boundedAlpha2 = min(max(alpha2, scalar(0)), scalar(1));
boundedAlpha3 = min(max(alpha3, scalar(0)), scalar(1));

if (phaseChange)
{

Info << "Phase Change \n" << endl;

dmdt = ((T - meltingT)/Lf)*mixture.cp();//*(1.0 - alpha1);

    forAll(dmdt, cellI)
    {
        if(alphaGeo)
        {
            if (!(alpha1[cellI] > 0.25) && ((alpha2[cellI] >= 0.5) || (alpha3[cellI] >= 0.5)))
                dmdt[cellI] *= 1;
            else
    			dmdt[cellI] *= 0;
        }      
        else
        {
            if ((alpha2[cellI] >= 0.25) || (alpha3[cellI] >= 0.25))
                dmdt[cellI] *= 1;
            else
                dmdt[cellI] *= 0;
        }
        Su[cellI] = dmdt[cellI]*dmdtRelaxAlpha.value();
        dmdt[cellI] = -dmdt[cellI]*dmdtRelaxAlpha.value();
    }
    
    dmdtSp = 0.0;
    zeroField Sp;
}
else
{
    dmdt = 0.0;
    dmdtSp = 0.0;
    zeroField Su;
    zeroField Sp;
}

zeroField divU;
