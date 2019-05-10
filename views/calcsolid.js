// Globale Variablen
var strucType = "A0";
var maxCompund = 0; 
var numAtom = 0;
var noAtom = 0;
var numPos = 0;
var lmsum = 0;
var crystForm = "gas";
var atomCharge = [];
var atomPos = [];
var espdofz = [];
var e_at = [];
var nvalence = 0;
//
// for plotting
var blable = ["1","2","3","4"]

//
// global physical paramters
var kboltz = 8.617e-5   // Boltzmannkonstante in eV/K
var rgas = 8.314 // universelle gaskonstante in J/(mol*K)
var ehart = 0.03675 // 1/Hartree-Energie
var anull = 0.529 // Bohr-Radius
var pi = 3.1415926;
//
// Unit transformation constants
var gpa = 160.2  // Umrechnung ev/Angström^3 auf GPa
var kgprol = 1.6606 
var mpros = 1000
var evnm = 1240
var evinkjmol = 96.485
var blin = 135 // emp. Faktor für Lindemannsche Schmelzpktformel


function setStructType(getType,getNum) {
    var wrtline = document.createTextNode("selected " + getType + " with " + getNum+ " Atoms.");
    var element = document.getElementById('StructureType').appendChild(wrtline);  
    
    maxCompund = getNum;
    strucType = getType;
    
    if (strucType == "Ah") {
        numAtom = 1;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]};           
    }


    if (strucType == "A1") {
        numAtom = 1;
        crystForm = "fcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]};           
    }


    if (strucType == "A3") {
        numAtom = 2;
        crystForm = "hex";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.000 , 0.000, 0.000]};      
        atomPos[2] =  { "Sort": 1 , "Pos":[ 0.667 , 0.333, 0.500]};               
    }

    if (strucType == "A2") {
        numAtom = 1;
        crystForm = "bcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]};           
    }    
    

    if (strucType == "A4") {
        numAtom = 2;
        crystForm = "fcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]}; 
        atomPos[2] =  { "Sort": 1 , "Pos":[ 0.25 , 0.25, 0.25]};           
    }
    
    if (strucType == "B1") {
        numAtom = 2;
        crystForm = "fcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]}; 
        atomPos[2] =  { "Sort": 2 , "Pos":[ 0.5 , 0.5, 0.5]};           
    }

    if (strucType == "B2") {
        numAtom = 2;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]}; 
        atomPos[2] =  { "Sort": 2 , "Pos":[ 0.5 , 0.5, 0.5]};           
    }

    if (strucType == "B3") {
        numAtom = 2;
        crystForm = "fcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.00 , 0.00, 0.00]}; 
        atomPos[2] =  { "Sort": 2 , "Pos":[ 0.25 , 0.25, 0.25]};           
    }


    if (strucType == "Gas") {
        numAtom = 2;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.50 , 0.50, 0.50]};  
        atomPos[2] =  { "Sort": 2 , "Pos":[ 0.51 , 0.51, 0.51]};             
    }       

    if (numAtom == 0) {
        var wrtline = document.createTextNode("selected " + getType + " still not configured :-(");
        var element = document.getElementById('StructureType').appendChild(wrtline);  
        var wrtline = document.createTextNode("use default simple cubic");
        var element = document.getElementById('StructureType').appendChild(wrtline); 
        numAtom = 1;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]};   
    }

    if (numAtom >= 1) {
        document.getElementById('selectAtoms').disabled = false; 
    }

}



function buildCrystal() {
     var rt = 1; 
     var eatom = 0.0; etotatom = 0.0; nvalat = 0; 

    var elemstr = "ReportCalc"
    var txtstr = "We calculate a crystal of structure type " + strucType 
                + "with " + numAtom + " atom(s) in the unit cell. " 
                + "The atoms are of " + maxCompund + " different compund(s)."
                + " The crystal form is " + crystForm;
    txtout(txtstr, elemstr);

    console.log(maxCompund);
     for (iatom = 1; iatom <= numAtom; iatom++ ) {
        e_at[iatom] = {"No":1, "sort:": 1, "erg": [0,0,0], "lmax": 2}; 
        console.log (espdofz[iatom]);

        isort = atomPos[iatom].Sort
        atomSort = espdofz[isort].atName;
        e_at[iatom].erg[0] = espdofz[isort].e_spd[0];
        e_at[iatom].erg[1] = espdofz[isort].e_spd[1];
        e_at[iatom].erg[2] = espdofz[isort].e_spd[2];              
        eatom =   e_at[iatom].erg[0]*espdofz[isort].q_spd[0] 
                + e_at[iatom].erg[1]*espdofz[isort].q_spd[1] 
                + e_at[iatom].erg[2]*espdofz[isort].q_spd[2];
        nvalat = espdofz[isort].q_spd[0] + espdofz[isort].q_spd[1] +espdofz[isort].q_spd[2];   
        etotatom = etotatom + eatom;
        nvalence = nvalence + nvalat; 
        if (e_at[iatom].erg[2] < 0) {
            e_at[iatom].lmax  = 3;
        }
        lmsum = lmsum +  power(e_at[iatom].lmax,2);  

        var elemstr = "ReportCalc";
        var txtstr = "The crystal has a " + atomSort 
                        + " Atom with valence energy " 
                        + form(eatom,3) + " ryd "
                        +"using lmax " + e_at[iatom].lmax + ".";
        txtout(txtstr,elemstr);
    }
    etotatom = form(etotatom,3);

    var elemstr = "ReportCalc";
    var txtstr = "The total atomic valence energy is " + etotatom 
                    + " ryd. We expecting " + lmsum 
                    + " eigenstates per k.";
    txtout(txtstr,elemstr);


     console.log (atomPos);
     console.log (atomPos[1].Pos);
     console.log (atomPos[1].Pos[1]);
     console.log (espdofz);
     console.log (e_at);
     console.log (e_at[1].erg);

     setCrystalValues();
     

     return rt;
}

function setCrystalValues() {

var GPoint = [0.0, 0.0, 0.0];

if (crystForm == "cub") {
   var t1 = [1.0, 0.0, 0.0];
   var t2 = [0.0, 1.0, 0.0];
   var t3 = [0.0, 0.0, 1.0];
   var APoint = [0.5, 0.5, 0.0];  // M
   var BPoint = [0.5, 0.0, 0.0];  // X
   var CPoint = [0.5, 0.5, 0.5];  // R
   var blable = ["M","G","X","R"];   
   var b1titel = "Bandindex M - G - X - R";
}

if (crystForm == "fcc") {
    var t1 = [0.0, 0.5, 0.5];
    var t2 = [0.5, 0.0, 0.5];
    var t3 = [0.5, 0.5, 0.0];
    var APoint = [0.5, 0.5, 0.5];  // L
    var BPoint = [0.0, 1.0, 0.0];  // X
    var CPoint = [0.5, 0.5, 0.0];  // K
    var blable = ["L","G","X","K"];
    var b1titel = "Bandindex L - G - X - K";
 }

 if (crystForm == "bcc") {
    var t1 = [-0.5, 0.5, 0.5];
    var t2 = [0.5, -0.5, 0.5];
    var t3 = [0.5, 0.5, -0.5];
    var APoint = [0.25, 0.25, 0.25];  // P
    var BPoint = [0.00, 0.50, 0.00];  // H
    var CPoint = [0.25, 0.25, 0.00];  // N
    var blable = ["P","G","H","N"];
    var b1titel = "Bandindex P - G - H - N";
 }

 var elemstr = "BandCalc";
 var txtstr = "The bandstructure is determined for the symmetry points " 
                 + "[" + APoint + "], "
                 + "[" + BPoint + "], "
                 + "[" + CPoint + "]"                                 
                 + " and of course the Gamma-Point G. "
                 + "The k-points are labeld as " + b1titel + ".";
 txtout(txtstr,elemstr);


//
// determine the real spacelattice vectors
var T = [t1,t2,t3];
var volume = numeric.det(T);

var elemstr = "ReportCalc";
var txtstr = "The real lattice vectors T are defined as " 
                + "{ [" + t1 + "] / [" + t2 + "] / [" + t3 + "] }.";
txtout(txtstr,elemstr);


//
// Determine the set of reciprocal space and k-Vectors 
var B = numeric.inv(T)

var elemstr = "ReportCalc";
var txtstr = "Which leads to following reciprocal space vectors B = " 
                + "{ [" + B + "]}.";
txtout(txtstr,elemstr);

var nkband = 4; var kpoint = [APoint , GPoint , BPoint , CPoint];
for (irun=0; irun < nkband; irun++) {
//    kpoint[irun] = [sum(kpoint(i,1).*b1),sum(kpoint(i,2).*b2),sum(kpoint(i,3).*b3)]
   kpoint[irun] = [numeric.sum(numeric.dot(kpoint[irun][0],B[0])),
                    numeric.sum(numeric.dot(kpoint[irun][1],B[1])),
                    numeric.sum(numeric.dot(kpoint[irun][2],B[2])),
                    ];

}

//
// build a super cell and determine the vectors between the neighbored cells
var ncut = 1; iq = 0; iat = 0; jat = 0;
var vecpos = []; vp=[];
for (h=-ncut; h <= ncut; h++) { 
    for (k=-ncut; k <= ncut; k++) { 
        for (l=-ncut; l <= ncut; l++) {
            iq = iq + 1
            vp = [h,k,l];
            vecpos[iq] = numeric.dot(T,vp);       
         } 
    }
 } 
 var nq = iq 
//
// determine the real atomic positions within the unit cell and the super cell
var xatomcor = [];
var atomcor = [];  
for (iat=1; iat <= numAtom; iat++) { 
    atomcor[iat] = numeric.dot(atomPos[iat].Pos,T);
}

console.log ("AtomCor ",atomcor);
igo = 0;
 for (iat=1; iat <= numAtom; iat++) {
    for (iq = 1; iq <= nq; iq++) {
       igo = igo + 1;
       vp = numeric.add(vecpos[iq],atomcor[iat]); 
       xatomcor[igo] =  { "xSort": iat , "qVec": iq, "xPos": vp};       
    }
}
numPos = igo;
//
// Determine the distances between the atoms
var dmin = numeric.rep([numAtom,numAtom],10000); totmin = 10000;
var near = []; igo = 0;

for (iat = 1; iat <= numAtom; iat++) {
    for (igo = 1; igo <= numPos; igo++) {
          jat = xatomcor[igo].xSort; 
          vq = numeric.sub(atomcor[iat],xatomcor[igo].xPos);            
          dist = sqrt(vq[0]*vq[0]+vq[1]*vq[1]+vq[2]*vq[2]);
            if (dist > 0) {
                if (dmin[iat-1][jat-1] > dist) {
                    dmin[iat-1][jat-1] = dist;
                    if (totmin > dmin[iat-1][jat-1]) {
                        totmin = dmin[iat-1][jat-1];                        
                    }
                }
            }
        }
}  
console.log("NN-Distance",totmin);

//   
// ...and the nearest neighbours
var sr = []; srd = 0; iaa = 0
for (iat = 1; iat <= numAtom; iat++) {
    near[iat] = 0;
    for (igo = 1; igo <= numPos; igo++) {
        iaa = iaa + 1;
        vq = numeric.sub(atomcor[iat],xatomcor[igo].xPos);            
        dist = sqrt(vq[0]*vq[0]+vq[1]*vq[1]+vq[2]*vq[2]);
        srd = 0;
        if (dist > 0) {     
                srd = power(totmin/dist,12);
                near[iat] = near[iat] + srd;    
            }   
        sr[iaa] = {"Atm1": iat, 
                     "Atm2":xatomcor[igo].xSort,
                     "qVec":xatomcor[igo].qVec,
                     "srd": srd, 
                     "dd": dist,
                     "dx": vq[0],
                     "dy": vq[1],
                     "dz": vq[2]};            
    }
}

//
// No we determine the lattice constant by Muffin-Tin Approx

console.log("minima",dmin);
var aest = numeric.rep([numAtom,numAtom],0.0);
for (iat = 1; iat <= numAtom; iat++) {
    isort = atomPos[iat].Sort
    for (jat = 1; jat <= numAtom; jat++) {
        jsort = atomPos[jat].Sort
         aest[iat-1][jat-1] = (espdofz[isort].ratom +espdofz[jsort].ratom)/dmin[iat-1][jat-1];
    }
}
console.log("aest",aest);
alat = maxM(aest)/100/anull;
console.log("Lattice Constant:",alat*anull);

var elemstr = "ReportCalc";
var txtstr = "From the investigation of the super cell of 27 neighbour unit cells, " 
            + "we determine that the smallest distance is " + form(totmin,3) + " fractional units. "
            + "This leads to an effective coordination number of  " + form(near[1],1) 
            + " for the first atom etc. Under the assumptions, that the neighbour atoms "
            + "are touchuing each other according to their atomic radius we can determine "
            + "the lattice constant a = " + form(alat*anull*100,1) + " pm."; 
txtout(txtstr,elemstr);

//
// now we define the tight biding parameters
var pdr= numeric.identity(3);  // power of distance dependenice (1/r^pdr)
var Vsig = numeric.rep([3,3],0) // sigma-boonding overlap
var Vpi0= numeric.rep([3,3],0) // pi-bonding overlap
var Ylm = numeric.rep([3,5],0) // speherical harmonics
var Elmg00r = numeric.rep([lmsum,lmsum],0) // part of matrix elemnetsvar 
var Elmg00i = numeric.rep([lmsum,lmsum],0) // part of matrix elemnets
var Elmg00 = numeric.rep([lmsum,lmsum],0) // part of matrix elemnets

var eev = numeric.rep([nkband,lmsum],0); Solv = [];
var psk= numeric.rep([9,9],0)  // slater-koster overlaps

console.log(Elmg00,numAtom,lmsum);
console.log(pdr);
console.log(Vsig,Vpi0,Ylm);
//
// DIstance dependencie
pdr[0][0] = 2; pdr[1][1] = 2; pdr[0][1] = 2; pdr[1][0] = 2;
pdr[2][0] = 3.5; pdr[2][1] = 3.5; pdr[0][2] = 3.5; pdr[1][2] = 3.5; 
pdr[2][2] = 5;    
//
// Overlapping parameters
var nss = -1.32; nsp = 1.42; nxx =  2.22; nxy = -0.63 
var ndds = -35.6; nddp = 19.2; nsd = -4.68; npds = -4.37; npdp = 2.02

Vsig[0][0] =  nss/power(alat,pdr[0][0]);
Vsig[0][1] =  nsp/power(alat,pdr[0][1]);
Vsig[1][0] =  -nsp/power(alat,pdr[1][0]);
Vsig[1][1] =  nxx/power(alat,pdr[1][1]);
Vsig[0][2] =  nsd/power(alat,pdr[0][2]);
Vsig[2][0] =  -nsd/power(alat,pdr[2][0]); 
Vsig[2][1] =  npds/power(alat,pdr[2][1]);
Vsig[1][2] =  -npds/power(alat,pdr[1][2]);
Vsig[2][2] =  ndds/power(alat,pdr[2][2]);
Vpi0[1][1] =  nxy/power(alat,pdr[1][1]);
Vpi0[1][2] =  npdp/power(alat,pdr[1][2]);
Vpi0[2][1] =  npdp/power(alat,pdr[2][1]);
Vpi0[2][2] =  nddp/power(alat,pdr[2][2]);

console.log("Vsig,Vpi",Vsig,Vpi0);

var kv = []; g00 = []
for (ik=0; ik < nkband; ik++) {
    iaa = 0;
    kv = kpoint[ik];
    console.log ("k-values",ik, kpoint[ik]);
    var kxp = 2*pi*kv[0]; kyp = 2*pi*kv[1]; kzp = 2*pi*kv[2];
    console.log ("k-values",ik,kxp,kyp,kzp);
//
// set the one-electrone nergies on the diagonals
    Elmg00r = numeric.rep([lmsum,lmsum],0)
    Elmg00i = numeric.rep([lmsum,lmsum],0)    
    var ialm = 0;
    for (iat = 1; iat <= numAtom; iat++) {
        for (il = 1; il <= e_at[iat].lmax; il++ ) {
            var ilm = 2*il - 1;
            for (var im = 1; im <= ilm; im++) {
                Elmg00r[ialm][ialm] = e_at[iat].erg[il-1];
                console.log("state",ialm,Elmg00r[ialm,ialm]);
                ialm = ialm +1;
            }
        }
    }    

    for (iat = 1; iat <= numAtom; iat++) {
        for (igo = 1; igo <= numPos; igo++) { 
        iaa = iaa + 1; 
           if (sr[iaa].srd > 0.01) {            
                if (sr[iaa].dd > 0) {
                    g00 = expi((sr[iaa].dx*kxp+sr[iaa].dy*kyp+sr[iaa].dz*kzp));                    
                    dl = sr[iaa].dx/sr[iaa].dd;
                    dm = sr[iaa].dy/sr[iaa].dd;
                    dn = sr[iaa].dz/sr[iaa].dd;

                    Ylm[0][0] = 1                        // s
                    Ylm[1][0] = dl                       // px
                    Ylm[1][1] = dm                       // py
                    Ylm[1][2] = dn                       // pz
                    Ylm[2][0] = sqrt(3)/2*(dl*dl-dm*dm)  // dx2-y2
                    Ylm[2][1] = dn*dn-(dm*dm+dl*dl)/2    // d3z2-r2
                    Ylm[2][2] = sqrt(3)*dl*dm            // dxy
                    Ylm[2][3] = sqrt(3)*dl*dn            // dxz
                    Ylm[2][4] = sqrt(3)*dm*dn            // dzy    

                    var imat = 0;
                    var ilv = 0;  
                    console.log ("bloch",g00,dl,dm,dn,sr[iaa].dd);                       
                    console.log ("harmonics",Ylm)                    
//                   console.log("Atom:",e_at);                
                    for (var il=0; il < e_at[iat].lmax; il++ ) {
                        var ilm = 2*(il+1) - 1;
                        for (var im = 0; im < ilm; im++) {
                            imat = imat + 1;
                            ilv = ilv + 1;
                            jat = xatomcor[igo].xSort;
                            var jmat = 0;
                            var jlv = 0;                
                            for (var jl=0; jl < e_at[jat].lmax; jl++) {
                                jlm = 2*(jl+1) - 1;
                                for (jm=0; jm < jlm; jm++) {
                                    jmat = jmat + 1;
                                    jlv = jlv + 1;
                                    pv = isequal(il,jl)*isequal(im,jm);
                                    dr = 1/power(sr[iaa].dd,pdr[il][jl]);
                                    ialm = (iat-1)*power(e_at[iat].lmax,2) + imat - 1;
                                    jalm = (jat-1)*power(e_at[jat].lmax,2) + jmat - 1;
                                    Esig = dr*sr[iaa].srd*Ylm[il][im]*Ylm[jl][jm]*Vsig[il][jl]*g00[1]*sqrt((4/near[iat]));
                                    Epi0 = dr*sr[iaa].srd*(pv - Ylm[il][im]*Ylm[jl][jm])*Vpi0[il][jl]*g00[1]; 
                                    Elmg00r[ialm][jalm] = Elmg00r[ialm][jalm] + Esig + Epi0;
                                    Esigi = dr*sr[iaa].srd*Ylm[il][im]*Ylm[jl][jm]*Vsig[il][jl]*g00[2]*sqrt((4/near[iat]));
                                    Epi0i = dr*sr[iaa].srd*(pv - Ylm[il][im]*Ylm[jl][jm])*Vpi0[il][jl]*g00[2]; 
                                    Elmg00i[ialm][jalm] = Elmg00i[ialm][jalm] + Esigi + Epi0i;                             
                             //       console.log (imat,jmat,Esig,Epi0);
                                    if (il == 1) { if (im == 2) { if (jl == 1){ if (jm == 0) { 
                    
                                    console.log ("Epi",Epi0,"dr",dr,"sr",sr[iaa].srd,"py1y2",(pv - Ylm[il][im]*Ylm[jl][jm]),"Vs",Vpi0[il][jl],"g00",g00[1],"sqrn",sqrt((4/near[iat])));                                    
                                    }}}}
                                    //                                        Elmg00(iat,ilv,jat,jlv,ik) = Elmg00(iat,ilv,jat,jlv,ik) + dr*sr(iat,jat,iq)*psk(ilv,jlv)*Vpi0(il,jl)*g00 //*(4/near(iat))^0.5                                      

//console.log(imat,jmat,iat,jat,ialm,jalm,Esig,Epi0);  
//console.log("dr",dr,"sr-dd",sr[iaa].dd,"ylm",Ylm[il][im]*Ylm[jl][jm],"Vsig",Vsig[il][jl],"Vpi",Vpi0[il][jl],"g00",g00[1],"damp",sqrt((4/near[iat])));                                                                      
                                }
                            }
                        }
                    }                      


                }
            }
        }
    }
    console.log ("Matrixelement Rea",Elmg00r);
    console.log ("Matrixelement Img",Elmg00i); 
    console.log ("Matrixelement tot",Elmg00);    
    var sumup = 0; sumdn = 0; 
    var Elmchk = numeric.sub(Elmg00r,Elmg00i);
    for (igo = 0; igo < lmsum; igo++) {
            for (jgo = igo+1; jgo < lmsum; jgo++ ) {
                sumup = sumup + Math.abs(Elmchk[igo][jgo]);
                sumdn = sumdn + Math.abs(Elmchk[jgo][igo]);
            }
    }
    console.log ("test Upper",sumup,sumdn,Elmchk)
    var corr = 0.0;
    if (sumup*sumdn < 0.001 ) {
        corr = sumaM(Elmg00i)/lmsum/lmsum*2;
        Elmg00i = numeric.rep([lmsum,lmsum],0);
    }


    Elmg00 = numeric.add(Elmg00r,Elmg00i);
    var Solv = numeric.eig(Elmg00); 
    var eeigr = Solv.lambda.x;
    var eeigi = Solv.lambda.y;
    if (eeigi == null) {
        eeigi = numeric.rep([lmsum,],0);
    }
    console.log ("Matrixelement er,ei",eeigr,eeigi);    

    var eeig = sortV(numeric.add(eeigr,eeigi));     
    for (igo = 0; igo < lmsum; igo++) {
            eeig[igo] = eeig[igo] - corr;
            if (igo > lmsum/2-1) {
                eeig[igo] = eeig[igo] + 2*corr;
            }
    }
    console.log ("Matrixelement Erg",eeig,corr,Solv);
  

//    var Solvi = numeric.eig(Elmg00i); 
//    var eeigi = sortV(Solvi.lambda.x);
//    console.log ("Matrixelement Ergi",eeigi,Solvi);
    
    for (irun = 0; irun < lmsum; irun++) {
        eev[ik][irun] = eeig[irun]/ehart;
    }
 //
 // transpose for plotting and analyzing   
    var eevk = numeric.transpose(eev);
    console.log ("Matrixelement eev",eev);
    console.log ("Matrixelement eevk",eevk);
//
// search for Fermi energy (rough) at Gamma-Point
    var nfermi = Math.floor(nvalence/2)-1;
    if (nfermi < 0) {
        nfermi = 0
    }
    if (isequal(nfermi,lmsum))
        nfermi = lmsum - 1
    }
    var gapdirect = eev[1][nfermi+1] - eev[1][nfermi];
    var efermi = (eev[1][nfermi+1] + eev[1][nfermi])/2;
    console.log ("Fermi-Energie",efermi,"nfermi",nfermi,"bandgap", gapdirect);

//
// Calculate the DOS by linear interpolation of the 4 k-Points
// G is the startpoint (in the middle of the BZ
// the other three pooints are on the plane on the BZ surface
var jdos=-1;
var ekkfit = []; ekf = []; efsort = [];
for (igo=0; igo < lmsum; igo++) {
    jdos = jdos + 1;    
    ekkfit[jdos] = eev[1][igo];
}
console.log("Egamma",ekkfit);
var nmix = 10000/lmsum/7;
for (irun =1; irun <= nmix; irun++) {
    kpo = sqrt(irun/nmix);
    for (igo=0; igo < lmsum; igo++) {
        db1 = eev[0][igo] - eev[1][igo];
        db2 = eev[2][igo] - eev[1][igo];
        db3 = eev[3][igo] - eev[1][igo];
        ekf[1] = eev[1][igo] +db1*kpo; 
        ekf[2] = eev[1][igo] +db2*kpo; 
        ekf[3] = eev[1][igo] +db3*kpo; 
        ekf[4] = (ekf[1] + ekf[2] + ekf[3])/3;
        ekf[5] = (ekf[2] + ekf[3])/2;
        ekf[6] = (ekf[1] + ekf[3])/2;
        ekf[7] = (ekf[1] + ekf[2])/2;        
        for (ia=1; ia <=7; ia++) {
            jdos = jdos + 1;
            ekkfit[jdos] = ekf[ia];
        }
    }
}

//
// Fermi again but now more accurate
var ndos = jdos;
var efsort = sortV(ekkfit);
var nfermi = Math.floor(ndos*nvalence/lmsum/2)
if (isequal(nfermi,0)) {
   nfermi = 1
}
if (isequal(nfermi,ndos)) {
   nfermi = ndos - 1
}
var efermi = (efsort[nfermi] + efsort[nfermi+1])/2;
var bandgap = efsort[nfermi+1]-efsort[nfermi];
var bandwidth = efsort[1] - efermi + bandgap/2
console.log("DOS", ekkfit);
console.log("EFermi", efermi,"bandgap",bandgap,"Bandwidth", bandwidth);
//
// Set the Fermi-Energy to zero
efsort = numeric.sub(efsort,efermi);
eevk = numeric.sub(eevk,efermi);
console.log("DOS", efsort);
var emaxx = maxV(efsort);
var eminx = minV(efsort);
var eminy = Math.floor(eminx/5)*5;
var emaxy = Math.ceil(emaxx/5)*5;
//
// Calulate the DOS withni a mesh of approriate size
console.log("Emins", eminx,eminy,"Emaxs",emaxx,emaxy);
var mesh = 100;
var defbox = (emaxx-eminx)/mesh;
var dosfit = hist(efsort,mesh);
normdos = sumV(dosfit)/(lmsum*2)*defbox;
dosfit = numeric.div(dosfit,normdos);
console.log ("Norm:",normdos,sumV(dosfit),sumV(dosfit)*defbox);
//
//  For plotting set energy labels 
var dlabel = []; dlabel[0] = form(eminy,1);
for (igo = 0; igo < mesh; igo++) {
    dlabel[igo] = ""; 
}
var nlabel = Math.floor((emaxy-eminy)/5);
var ntebox = Math.floor(mesh/nlabel);
for (igo = 1; igo <= nlabel; igo++) {
    dlabel[igo*ntebox] = Math.floor(eminy+5*igo);
}
dlabel[0] = form(eminy,1);
//dlabel[mesh-1] = form(emaxy,1);
console.log ("DOS",dosfit,dlabel);

//    var  A = [[1,0,1],[0,2,0],[0.3,0,1]]; Am = A;
//    var u  = sqrt(numeric.norm2(A));
//    var rn = numeric.diag(numeric.getDiag(A));
//    console.log (A,u, rn);
//    for (irun = 1; irun < 5; irun++) {
//        Am = numeric.dot(rn,A);
//        u  = sqrt(numeric.norm2(Am));
//       rn = numeric.getDiag(numeric.div(Am,u));
//       console.log (A, Am, u, rn);
//       rn = numeric.diag(rn);
//   }

// plot Bandstructure
var lineChartData = {
    labels: eval(blable),
    datasets: [
        {
            fillColor: "rgba(255, 255, 255, 0)",
            strokeColor: "rgba(000,000,000,1)",
            data: eval(eevk[0])
        }
    ]
}

for (il = 1; il <= lmsum; il++) {
    lineChartData.datasets[il] = {"data": eevk[il] };
}
var myLine = new Chart(document.getElementById("bandstructure").getContext("2d")).Line(lineChartData);

var barChartData = {
    labels : eval(dlabel),
    datasets : [
        {
            fillColor : "rgba(000,000,000,1)",
            strokeColor : "rgba(45,53,65,1)",
//            highlightFill: "rgba(220,220,220,0.75)",
//            highlightStroke: "rgba(220,220,220,1)",
            data : eval(dosfit)
        }
    ]

}

var myLine = new Chart(document.getElementById("densityofstates").getContext("2d")).Bar(barChartData);





}