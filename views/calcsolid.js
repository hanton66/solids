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
var nexpanh = 0;
var compoundname = " ";
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


function setStructType(getType) {

    strucType = getType;
    
    if (strucType == "Ah") {
        maxCompund = 1;
        numAtom = 1;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]};           
    }


    if (strucType == "A1") {
        maxCompund = 1;
        numAtom = 1;
        crystForm = "fcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]};           
    }

    if (strucType == "A2") {
        maxCompund = 1;
        numAtom = 1;
        crystForm = "bcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]};           
    }    
    

    if (strucType == "A3") {
        maxCompund = 1;
        numAtom = 2;
        crystForm = "hex";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.000 , 0.000, 0.000]};      
        atomPos[2] =  { "Sort": 1 , "Pos":[ 0.667 , 0.333, 0.500]};               
    }


    if (strucType == "A4") {
        maxCompund = 1;
        numAtom = 2;
        crystForm = "fcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]}; 
        atomPos[2] =  { "Sort": 1 , "Pos":[ 0.25 , 0.25, 0.25]};           
    }

    if (strucType == "A9") {
        maxCompund = 1;
        numAtom = 2;
        crystForm = "hex";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.000 , 0.000, 0.000]};   
        atomPos[2] =  { "Sort": 1 , "Pos":[ 0.333 , 0.667, 0.005]};                   
    }
    
    if (strucType == "B1") {
        maxCompund = 2;
        numAtom = 2;
        crystForm = "fcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]}; 
        atomPos[2] =  { "Sort": 2 , "Pos":[ 0.5 , 0.5, 0.5]};           
    }

    if (strucType == "B2") {
        maxCompund = 2;
        numAtom = 2;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]}; 
        atomPos[2] =  { "Sort": 2 , "Pos":[ 0.5 , 0.5, 0.5]};           
    }

    if (strucType == "B3") {
        maxCompund = 2;
        numAtom = 2;
        crystForm = "fcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.00 , 0.00, 0.00]}; 
        atomPos[2] =  { "Sort": 2 , "Pos":[ 0.25 , 0.25, 0.25]};           
    }

    if (strucType == "C3") {
        maxCompund = 2;
        numAtom = 6;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.25 , 0.25, 0.25]}; 
        atomPos[2] =  { "Sort": 1 , "Pos":[ 0.25 , 0.75, 0.75]}; 
        atomPos[3] =  { "Sort": 1 , "Pos":[ 0.75 , 0.25, 0.75]}; 
        atomPos[4] =  { "Sort": 1 , "Pos":[ 0.75 , 0.75, 0.25]};
        atomPos[5] =  { "Sort": 2 , "Pos":[ 0.00 , 0.00, 0.00]}; 
        atomPos[6] =  { "Sort": 2 , "Pos":[ 0.50 , 0.50, 0.50]};                  
    }

    if (strucType == "C9") {
        maxCompund = 2;
        numAtom = 6;
        crystForm = "fcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.000 , 0.000, 0.000]}; 
        atomPos[2] =  { "Sort": 1 , "Pos":[ 0.250 , 0.250, 0.250]}; 
        atomPos[3] =  { "Sort": 2 , "Pos":[ 0.125 , 0.125, 0.125]}; 
        atomPos[4] =  { "Sort": 2 , "Pos":[ 0.375 , 0.375, 0.125]};
        atomPos[5] =  { "Sort": 2 , "Pos":[ 0.375 , 0.125, 0.375]}; 
        atomPos[6] =  { "Sort": 2 , "Pos":[ 0.125 , 0.375, 0.375]};                  
    }

    if (strucType == "E21") {
        maxCompund = 3;
        numAtom = 5;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.000 , 0.000, 0.000]}; 
        atomPos[2] =  { "Sort": 2 , "Pos":[ 0.500 , 0.500, 0.500]}; 
        atomPos[3] =  { "Sort": 3 , "Pos":[ 0.000 , 0.500, 0.500]}; 
        atomPos[4] =  { "Sort": 3 , "Pos":[ 0.500 , 0.000, 0.500]};   
        atomPos[5] =  { "Sort": 3 , "Pos":[ 0.500 , 0.500, 0.000]};                      
    }    


    if (strucType == "L12") {
        maxCompund = 2;
        numAtom = 4;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 2 , "Pos":[ 0.000 , 0.000, 0.000]}; 
        atomPos[2] =  { "Sort": 1 , "Pos":[ 0.000 , 0.500, 0.500]}; 
        atomPos[3] =  { "Sort": 1 , "Pos":[ 0.500 , 0.000, 0.500]}; 
        atomPos[4] =  { "Sort": 1 , "Pos":[ 0.500 , 0.500, 0.000]};                 
    }

    if (strucType == "V0") {
        maxCompund = 1;
        numAtom = 1;
        crystForm = "fcc";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.00 , 0.00, 0.00]};              
    } 

    if (strucType == "V1") {
        maxCompund = 1;
        numAtom = 2;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.49 , 0.49, 0.49]};  
        atomPos[2] =  { "Sort": 1 , "Pos":[ 0.51 , 0.51, 0.51]};             
    } 

    if (strucType == "V2") {
        maxCompund = 2;
        numAtom = 2;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.48 , 0.48, 0.48]};  
        atomPos[2] =  { "Sort": 2 , "Pos":[ 0.52 , 0.52, 0.52]};             
    }       

    if (numAtom == 0) {
        var wrtline = document.createTextNode("selected " + getType + " still not configured :-(");
        var element = document.getElementById('StructureType').appendChild(wrtline);  
        var wrtline = document.createTextNode("use default simple cubic");
        var element = document.getElementById('StructureType').appendChild(wrtline); 
        numAtom = 1;
        maxCompund = 1;
        crystForm = "cub";
        atomPos[1] =  { "Sort": 1 , "Pos":[ 0.0 , 0.0, 0.0]};   
    }

    if (numAtom >= 1) {
        document.getElementById('selectAtom1').disabled = false; 
        document.getElementById('deleteAtom').disabled = false;         
    }

    if (numAtom == 0) {
        document.getElementById('selectAtom1').disabled = true; 
        document.getElementById('deleteAtom').disabled = true;         
    }   


    //var element = document.getElementById('StructureType');
   // element.parentNode.removeChild(element);
    var wrtline = document.createTextNode("selected " + getType + " with "  + maxCompund + " atoms." );
    var element = document.getElementById('StructureType');  
    element.replaceChild(wrtline, element.firstChild);      

}



function buildCrystal() {
     var rt = 1; 
     var eatom = 0.0; etotatom = 0.0; nvalat = 0; 

    var elemstr = "ReportCalc"
    var txtstr = " We calculate a crystal of structure type " + strucType 
                + " with " + numAtom + " atom(s) in the unit cell. " 
                + " The atoms are from " + maxCompund + " different compund(s). "
                + " The crystal form is " + crystForm;
    txtout(txtstr, elemstr);

    console.log(maxCompund);

    var compelem=[]; companza=[];

    var cwei = maxCompund/numAtom; cmin = 1000;
    for (igo=1; igo <= maxCompund; igo++) {
        companza[igo] = 0;
        for (iatom=1; iatom <= numAtom; iatom++) {
            isort = atomPos[iatom].Sort;
            if (isequal(igo,isort)) {
                compelem[isort] = espdofz[isort].atName;
                companza[isort] = companza[isort] + cwei;
            }
        }
        if (companza[igo] < cmin) {
            cmin = companza[igo];
        }
    }

    var cnum = "";
    for (igo=1; igo <= maxCompund; igo++) {
        cnum = form(companza[igo]/cmin,0);
        if (cnum == 1) { cnum = ""; }
        compoundname = compoundname + compelem[igo] + cnum;
    }

     for (iatom = 1; iatom <= numAtom; iatom++ ) {
        e_at[iatom] = {"No":1, "sort:": 1, "erg": [0,0,0], "lmax": 2, "pointer": lmsum}; 
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
        nexpanh  = nexpanh + espdofz[isort].nexp;
        if (e_at[iatom].erg[2] != 0) {
            e_at[iatom].lmax  = 3;
        }
        lmsum = lmsum +  power(e_at[iatom].lmax,2);  

        var elemstr = "ReportCalc";
        var txtstr = " The crystal has a " + atomSort 
                        + " Atom with a valence energy of  " 
                        + form(eatom,3) + " ryd "
                        +" using a maximum angular momentum of " + e_at[iatom].lmax + " ."
                        +" The fractional atom position in the unit cell is at { " 
                        + atomPos[iatom].Pos + " } .";
        txtout(txtstr,elemstr);
    }
    etotatom = form(etotatom,3);

    var elemstr = "ReportCalc";
    var txtstr = "The total atomic valence energy is " + etotatom 
                    + " ryd. We expecting " + lmsum 
                    + " eigenstates per k-vector. ";
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

 if (crystForm == "hex") {
    var cbya = sqrt(8/3);
    var t1 = [0.5, form(-sqrt(3)/2,3), 0.0];
    var t2 = [0.5,  form(sqrt(3)/2,3), 0.0];
    var t3 = [0.0, 0.0, form(cbya,3)*2];
    var APoint = [0.50, 0.50, 0.00];  // K
    var BPoint = [0.00, 0.00, 0.50];  // A
    var CPoint = [0.50, 0.00, 0.50];  // L
    var blable = ["K","G","A","L"];
    var b1titel = "Bandindex K - G - A - L";
 }


 var elemstr = "BandCalc";
 var txtstr = " The bandstructure is determined for the following symmetry points " 
                 + "[" + APoint + "], "
                 + "[" + BPoint + "], "
                 + "[" + CPoint + "]"                                 
                 + " and the Gamma-Point G. "
                 +  "The k-points are labeld as " + b1titel + " .";
 txtout(txtstr,elemstr);


//
// determine the real spacelattice vectors
var T = [t1,t2,t3];
var volume = numeric.det(T);

var elemstr = "ReportCalc";
txtout("",elemstr);
var txtstr = " The real space lattice vectors R are defined as " 
                + "{ [" + t1 + "] / [" + t2 + "] / [" + t3 + "] }.";
txtout(txtstr,elemstr);


//
// Determine the set of reciprocal space and k-Vectors 
var B = numeric.transpose(numeric.inv(T));
//var B = numeric.transpose(Z);
console.log ("B Matrix ",B);
console.log ("B-1 Vektor ",B[0]);

var elemstr = "ReportCalc";
var txtstr = " Which leads to the following reciprocal space vectors B = " 
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
console.log ("K-Points ",kpoint);

//
// build a super cell and determine the vectors between the neighbored cells
var ncut = 1; iq = 0; iat = 0; jat = 0;
var vecpos = []; vp=[];
for (h=-ncut; h <= ncut; h++) { 
    for (k=-ncut; k <= ncut; k++) { 
        for (l=-ncut; l <= ncut; l++) {
            iq = iq + 1
            vp = [h,k,l];
            vecpos[iq] = numeric.dot(vp,T);  
//            vecpos[iq] = [numeric.dot(h,T[0])),
//                          numeric.sum(numeric.dot(k,T[1])),
//                          numeric.sum(numeric.dot(l,T[2])),];     
         } 
    }
 } 
 var nq = iq 
 console.log ("Vectorpos ",vecpos);
//
// determine the real atomic positions within the unit cell and the super cell
var xatomcor = [];
var atomcor = [];  
for (iat=1; iat <= numAtom; iat++) { 
    atomcor[iat] = numeric.dot(atomPos[iat].Pos,T);
    if (isequal(crystForm,'fcc'))  { 
        atomcor[iat] = atomPos[iat].Pos;
     }
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
txtout("",elemstr);
var txtstr = " From the investigation of the super cell of 26 neighbour unit cells, " 
            + " we determine that the smallest distance is " + form(totmin,3) + " fractional units. "
            + " This leads to an effective coordination number of  " + form(near[1],0) 
            + " for the first atom etc. Under the assumption, that the neighbour atoms "
            + " are touching each other according to their atomic radius we can estimate "
            + " the lattice constant to a = " + form(alat*anull*100,1) + " pm."; 
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
//                    console.log ("bloch",g00,dl,dm,dn,sr[iaa].dd);                       
//                    console.log ("harmonics",Ylm)                    
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
                                    ialm = e_at[iat].pointer + imat - 1;
                                    jalm = e_at[jat].pointer + jmat - 1;
                                    Esig = dr*sr[iaa].srd*Ylm[il][im]*Ylm[jl][jm]*Vsig[il][jl]*g00[1]*sqrt((4/near[iat]));
                                    Epi0 = dr*sr[iaa].srd*(pv - Ylm[il][im]*Ylm[jl][jm])*Vpi0[il][jl]*g00[1]; 
                                    Elmg00r[ialm][jalm] = Elmg00r[ialm][jalm] + Esig + Epi0;
                                    Esigi = dr*sr[iaa].srd*Ylm[il][im]*Ylm[jl][jm]*Vsig[il][jl]*g00[2]*sqrt((4/near[iat]));
                                    Epi0i = dr*sr[iaa].srd*(pv - Ylm[il][im]*Ylm[jl][jm])*Vpi0[il][jl]*g00[2]; 
                                    Elmg00i[ialm][jalm] = Elmg00i[ialm][jalm] + Esigi + Epi0i;
//    in the case of Salte-Koster => dr*sr(iat,jat,iq)*psk(ilv,jlv)*Vpi0(il,jl)*g00                                      
                               }
                            }
                        }
                    }                      


                }
            }
        }
    }
//
// Now we solve the Eigen equation of the atomic states
// Unfortunately we can only calculate eigenvalues from real matirx.
// Therfore we gather the imaginery elements but them as shift to
// the valence band (-corr) and conduction band (+corr) for non-symm Matrices.
// For antisymmetric Matrix we assume that H* = real(H) + imag(H) delivers similar 
// eigenvalues like H
    var Elmg00 = numeric.add(Elmg00r,Elmg00i);
    var matSym = chkMatSy(Elmg00);
    if (matSym < 0.001) {
        var Solv = numeric.eig(Elmg00);
        var eeigr = Solv.lambda.x; 
        var eeigi = Solv.lambda.y;  
        if (eeigi == null) {
            eeigi = numeric.rep([lmsum],0);
        }           
        var eeig = sortV(numeric.add(eeigr,eeigi));                  
    }
    if (matSym > 0.001) {
        var Solv = numeric.eig(Elmg00r); 
        var eeigr = Solv.lambda.x;     
        var eeig = sortV(eeigr);                 
    }   

    var corr = 0.0;
    corr = sumaM(Elmg00i)/lmsum/lmsum;    
    for (igo = 0; igo < lmsum; igo++) {
            eeig[igo] = eeig[igo] - corr;
            if (igo > lmsum/2-1) {
                eeig[igo] = eeig[igo] + 3*corr;
            }
    }
    console.log ("Matrixelement Erg",eeig,corr,Solv);
    
    for (irun = 0; irun < lmsum; irun++) {
        eev[ik][irun] = eeig[irun]/ehart;
    }
 //
 // transpose energie values for plotting and analyzing   
    var eevk = numeric.transpose(eev);
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


    var elemstr = "BandCalc";
    txtout("",elemstr);
    var txtstr =  " The following line chart shows the bandstructure of the compound. "
                + " Occupied valence bands are in green. Unoccupied or not fully occupied " 
                + "conduction bands are plotted in red. "
                + " The energy scale in the plot is set relative to the Fermi-Energy (= 0 eV)"; 
    txtout(txtstr,elemstr);
    txtout("",elemstr);

    var elemstr = "BandCalc2";
    txtout("",elemstr);
    var txtstr =  " We evaluate the electronic structure of the Gamma-Point in the "
                + " middle of the Brillouin-Zone. The direct band gap is " + form(gapdirect,2)     
                + " eV. There are " + form(nfermi+1,0) + " valence bands occupied with " 
                + nvalence + " electrons. The Fermi energy is about " + form(efermi,2) + " eV. ";
    txtout(txtstr,elemstr);
    var txtstr =  " The deepest electronic band lies at " + form(eev[1][0]-efermi,2) + " eV."; 
    txtout(txtstr,elemstr);


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

var elemstr = "DosCalc";
txtout("",elemstr);
var txtstr =  " We approximate the density-of-state (DOS) for the compound by a linear interpolation "
            + " between the energy values of the calculated k-Point. We assume, that the four "
            + " k-points in our plot roughly describe the four corners of a representative Brillouin segment. "
            + " Withn this segment we run from G to the three remaining k-Points on the surface of the BZ-zone. "
            + " This is like going from the center of a sphere to the surface and meanns that the "
            + " weight of the calculated outer interpolated energy values are scaling with k-square, "
            + " whereas k means the absolute value of the k-vector. The following bar chart shows the "
            + " calculated DOS of the investigated compund.";
txtout(txtstr,elemstr);


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
// console.log("DOS", ekkfit);
// console.log("EFermi", efermi,"bandgap",bandgap,"Bandwidth", bandwidth);
//
// Set the Fermi-Energy to zero
efsort = numeric.sub(efsort,efermi);
eevk = numeric.sub(eevk,efermi);
// console.log("DOS", efsort);
var emaxx = maxV(efsort);
var eminx = minV(efsort);
var eminy = Math.floor(eminx/5)*5;
var emaxy = Math.ceil(emaxx/5)*5;
//
// Calulate the DOS withni a mesh of approriate size
// console.log("Emins", eminx,eminy,"Emaxs",emaxx,emaxy);
var mesh = 100;
//var defbox = (emaxx-eminx)/mesh;
var defbox = (emaxy-eminy)/mesh;
var dosfit = hist(efsort,mesh);
normdos = sumV(dosfit)/(lmsum*2)*defbox;
dosfit = numeric.div(dosfit,normdos);
console.log ("Norm:",normdos,sumV(dosfit),sumV(dosfit)*defbox);

//
// Calculate an energy mesh and labels for plotiing
//dlabel = getlabels(eminy,emaxy,mesh,5);
efspace = numeric.linspace(eminy,emaxy,mesh);
dlabel = getlabel(efspace,10);
//dlabel[mesh-1] = form(emaxy,1);
console.log ("DOS",dosfit,dlabel,efspace);

//
// Searching for the Fermi-Level in th DOS
var dosatef = 0;
for (igo = 0; igo < mesh-1; igo++ ) {
    if (efspace[igo] < 0) {
        if (efspace[igo+1] > 0) {
           dosatef = (dosfit[igo]+dosfit[igo+1])/2;
           nfermi = igo;
           console.log("DOS,nf",dosatef,nfermi)
        }
    }
    if (efspace[igo] == 0) {       
        dosatef = dosfit[igo];
        nfermi = igo;      
    }   
}

var gapratio = bandgap/gapdirect
var substtype = " metal";
var gaptype = "a";
if (isequal(dosatef,0)) {
    substtype = "n insulator";
    if (bandgap < 2.5) {
       substtype = " semiconductor";
        gaptype = "a direct";
        if (gapratio < 1) {
            gaptype = 'an indirect';
        }
    }
}


var elemstr = "DosCalc2";
txtout("",elemstr);
var txtstr =  " The DOS shows the calculated energy values in a range from " + form(eminx,2) 
            + " eV in the valence band to " + form(emaxx,2) + " eV in the conduction band.";
txtout(txtstr,elemstr);
var txtstr =  " We can now determine the Fermi-Energy more accurate because now all k-Points are "
            + " involved. The Fermi-Energy determined from the DOS is about " + form(efermi,2) + " eV. "
            + " The total bandwidth of the valence band can be determined to " + form(bandwidth,2) + " eV. "
            + " The band gap analysis of the DOS results in a gap of " + form(bandgap,2) + " eV."
            + " The charge density at the Fermi-Energy is " + form(dosatef,2) + " e-/eV. "
            + " The compund seems to be " + gaptype + substtype + ".";
txtout(txtstr,elemstr);


//
// Now we calculate the bonding energy in the crysatl according
// to the Fermi-Dirac Statistics (here at 1 Kelvin)
var eweight = [];
var evalence = 0; 
for (igo = 0; igo < mesh; igo++) {
    eweight[igo] = 1.0/(exp(efspace[igo]/(kboltz*1))+1);
    evalence = evalence + (efspace[igo]+efermi)*dosfit[igo]*eweight[igo];
}
evalence = evalence*defbox;
etotatom = etotatom/ehart;  // Changing from ryd to eV
ebond = evalence - etotatom;
console.log ("Valence bond enr",ebond )

//
// Now we calculate the Optical Properties from the DOS
// base of the optical properties is the (Pseudo-)Joint-DOS 
var pjdos=[]; hnuf=[];lambda=[];omega=[]; eps2=[]; eps1 = [];
var mesh = dosfit.length;
pjdos = numeric.rep([mesh],0);
for (jmesh = 0; jmesh < nfermi; jmesh++) {
    for (gmesh = nfermi; gmesh < mesh; gmesh++ ) {
        pjdos[gmesh-jmesh] = pjdos[gmesh-jmesh] + dosfit[jmesh]*dosfit[gmesh]*defbox*defbox;
    }
}
pjdos[0] = pjdos[1]/2;
console.log ("PJOS",pjdos);
//
// using a constant oscillator strength we calculated the transition probability eps2
// and via Kronig-Kramer relationship the dielectric/suceptibility functions
var alatA = alat*anull;
var volnull = volume*power(alatA,3);
var fosz = 22000/((lmsum - nvalence/2)*(nvalence/2));
var hnuf = numeric.linspace(defbox,emaxy-eminy+defbox,mesh);
console.log("ev-mesh",hnuf);
var olabel = getlabel(hnuf,10);
for (igo = 0; igo < mesh; igo++) {
    omega[igo] = form(hnuf[igo],2);
}
hmin = Math.round(0.8/defbox);
hmax = Math.round(12.4/defbox);
if (hmax > mesh) {hmax = mesh}
for (igo = 0; igo < mesh; igo++) {
    eps2[igo] = pi/volnull*fosz*pjdos[igo]/power(hnuf[igo],2)
}
var dshift = defbox/2;
for (igo = 0; igo < mesh; igo++) {
    omegas = omega[igo]+dshift;
    var sumopes = 0;
    for (jgo = 0; jgo < mesh; jgo++) {
        sumopes = sumopes + omega[jgo]*eps2[jgo]/(power(omega[jgo],2)-power(omegas,2))
    }
    eps1[igo] = 1 + 2/pi*sumopes;
}
//
// which give us the frquency dependend optical indices, reflection and transmission
nindex = sqrt(0.5*(sqrt(power(eps1[0],2)+power(eps2[0],2))+eps1[0]));
kindex = sqrt(0.5*(sqrt(power(eps1[0],2)+power(eps2[0],2))-eps1[0]));
reflection = power((nindex-1)/(nindex+1),2);
transmission = 1 - reflection;
//
// Norm eps1/2 for plotting
var normeps1 = []; normeps2=[];
for (igo = hmin; igo < hmax; igo++) {
    normeps2[igo-hmin] = eps2[igo]; 
    normeps1[igo-hmin] = eps1[igo]; 
    lambda[igo-hmin] = form(evnm/hnuf[igo],0);  
}
normeps1 = reverseV(scaleV(normeps1));
normeps2 = reverseV(scaleV(normeps2));
lambda = getlabel(reverseV(lambda),10);
//
// Calculate the optical gap using the Joint-DOS - could be very different in case of f-electrons
var optigap = defbox;
for (igo=1; igo <= mesh; igo++) {
    if (pjdos[mesh-igo] > 0.0001) {
        optigap = hnuf[mesh-igo];
//        console.log ("Optik Gap",optigap,defbox,igo,mesh-igo);
    }
}
var absorpedge = evnm/optigap;
if (isequal(optigap,defbox)) {
    optigap = 0;
}
var spektrumtyp = "extreme UV"
if (absorpedge > 200) {
    spektrumtyp = "near UV"; 
    if (absorpedge > 380) {
        spektrumtyp = "visible";
        if (absorpedge > 780) {
            spektrumtyp = "near IR";
            if (absorpedge > 3000) {
                spektrumtyp = "mid IR to microwave";
            }
        }
    }
}


var elemstr = "OptCalc";
txtout("",elemstr);
var txtstr =  " Now we try to approximate the optical properties of the compund. "
            + " Optical transitions occur from exciting electrons from the "
            + " valence bands to the conduction bands by light absorption. " 
            + " Therefore we can approximate the probaliity to excite an electron "
            + " by folding the valence band DOS to the conduction band DOS."
            + " E.g. we have a high transition probability, "
            + " if a photon provides exactly the energy difference between a DOS peak in the "
            + " valence band and a DOS peak in the conduction band. ";
txtout(txtstr,elemstr);
var txtstr =  " The folded DOS is called the Joint-Density-of-States. "
            + " In our approximation we do not take into account forbidden transitions "
            + " and calculate the JDOS directly from the DOS. Normally one has to make "
            + " this calculation for each k-Value seperately.";
txtout(txtstr,elemstr);


var elemstr = "OptCalc2";
txtout("",elemstr);
var txtstr =  " From the JDOS we can calculate the dielectric and suszeptibility functions by using the Kramer-Kronings relation [6]. "
            + " These are shown in the graph below where eps_1 is plotted in Magenta and eps2 is "
            + " plotted in Cyan. The eps2-function represents the absorption spectra of the compound. "
            + " The values of eps_1 and eps_2 are plotted against the wavelength of a photon in nm."
txtout(txtstr,elemstr);
var txtstr =  " From the eps1 and eps2 values at infinity we can calculate the refractive index [7]"
            + " which is n=" + form(nindex,2) + ", the refelction R= " + form(reflection*100,1)
            + " % and the transmission T= " + form(transmission,1)*100 + " %. ";
txtout(txtstr,elemstr);
var txtstr =  " Furthermore we can calculate the absorption edge " + form(absorpedge,1) +" nm"
            + " in eps2 and the optical band gap = " + form(optigap,2) + " eV."
            + " The abosrption edge lies in the " + spektrumtyp + " spectrum.";
txtout(txtstr,elemstr);


//console.log ("PJOS",pjdos,omega);
//console.log ("eps2",eps2,lambda,pi,volnull,fosz,hnuf);
//console.log ("Ep1",eps1)
//console.log ("n,k,R,T",nindex,kindex,reflection,transmission);
//console.log ("Plot e1,e2",lambda,normeps1,normeps2);
//console.log ("Opic",optigap,absorpedge);

//
// Calculate some physical Properties
var molmasse = 0;
for (iatom = 1; iatom <= numAtom; iatom++) {
    isort = atomPos[iatom].Sort
    var atomchg = espdofz[isort].atCharge;
    molmasse = molmasse + 0.0063*atomchg*atomchg + 2.0886*atomchg -  1.3453;
} 
var dichte = molmasse/volnull*kgprol;


var nexp = -nexpanh/numAtom;
var rwigner = power(3*volnull/4/pi/numAtom,1/3);
var bulkmodul = -gpa*ebond/volnull*power((nexp+1),2);
var betab = bulkmodul/gpa;
var shearmodul = 0.4*bulkmodul;
var c11modul = bulkmodul+4/3*shearmodul;
var youngmodul = 9*bulkmodul*shearmodul/(3*bulkmodul+shearmodul);
var debye = 177.5*power(numAtom,1/3)*power(volnull,1/6)*sqrt(bulkmodul/molmasse);
var vlwav = sqrt((3*bulkmodul+4*shearmodul)/(3*dichte))*mpros;
var vtwav = sqrt(shearmodul/dichte)*mpros;
var schall = power(1/3*(2/power(vtwav,3)+1/power(vlwav,3)),(-1/3));
var Enull = -betab*volnull/power(nexp+1,2)/numAtom;
var vaperg = -Enull*evinkjmol 
var rhoint=power(bulkmodul/1900,3/5);
var grueneisen= -1/6-nexp 
var hardness= 10*tanh(2*bulkmodul/500)+0.2
var vdiss =exp(-1/nexp+log(volnull))
var tschmelz = 3e-5*power(debye,2)*power(volnull/numAtom,2/3)*molmasse/numAtom
// var tschmelz = 3e-5*power(debye,2)*power((volnull/numAtom),(2/3))*molmasse/Nnu

var elemstr = "PhysCalc";
txtout("",elemstr);
var txtstr = " From the total valence energy calculation of the DOS we try to estimate a number of some"
            + " physical properties. First of all we can assume from the atomic composition same basic "
            + " static properties."; 
txtout(txtstr,elemstr);
txtout("",elemstr);
var txtstr = " The molare mass of the compound is about " 
                + form(molmasse,1) + " g/mol."
                + " The equlibrium volume of the unit cell is " + form(volnull,2) + " A^3. "
                + " That leads to a mass density of the compund of " + form(dichte,3) + " kg/l. ";
txtout(txtstr,elemstr);
txtout("",elemstr);
var txtstr = " Using an equation of state, we can determine from the energy/volume relation "
            + " other meachnical, thermal and elastic properties [4]. The bonding energy per unit cell "
            + " is about" + form(ebond,2) + " eV, which is roughly a dissoziation energy of "
            + form(Enull*evinkjmol,1) + " kj/mol.";
txtout(txtstr,elemstr);
var txtstr = " The bulk modulus is " + form(bulkmodul,1) + " GPa. The other eleastic constants [8]"
            + " can be approximated in the case of a polycrystalline substance. The shearmodul is "
            + form(shearmodul,1) + " GPa and the Young moduls is " + form(youngmodul,1) + " GPa."
            + " The tetragonal shear modul is " + form(c11modul,1) + " GPa.";
txtout(txtstr,elemstr);
var txtstr = " The Debye temperature [5] of this compound is " + form(debye,1) + " K. "
            + " The speed of sound is " + form(schall,1) + " m/s with a transversal velocity part "
            + " of " + form(vtwav,1) + " m/s and a logituidinal part of " + form(vlwav,1) + " m/s.";          ;
txtout(txtstr,elemstr);



console.log (rwigner, nexp);
console.log (bulkmodul, shearmodul, c11modul, youngmodul);
console.log (debye, vlwav,vtwav,schall,Enull, grueneisen, hardness,vdiss);

//
// now we calculate some thermal properties
var Tmax = 600; Tmin = 10; cvt=[]; xdt = []; mesh = 20;
var tspace = numeric.linspace(Tmin,Tmax,mesh);
// var tlabel = numeric.round(tspace);
var tlabel = getlabel(tspace,10);
var meanmol = molmasse/numAtom;
for (igo = 0; igo < mesh; igo++) {
    var sumxdt = 0;
    for (jgo = 0; jgo < mesh; jgo++) {
        xdt[jgo] = debye/tspace[igo]/1000 +  (debye/tspace[igo] - debye/tspace[igo]/1000)/(mesh-1)*jgo;
        sumxdt = sumxdt + power(xdt[jgo],4)*exp(xdt[jgo])/power(exp(xdt[jgo])-1,2);
    }
    dtx = xdt[1] - xdt[0]; 
    cvt[igo] = 9*rgas*power((tspace[igo]/debye),3)*sumxdt*dtx/meanmol*1000;
    if (tspace[igo] < 293) {
        var cvrt = cvt[igo];
    }
}
var alphal = grueneisen/3*cvrt/bulkmodul*dichte;

var elemstr = "PhysCalc2";
txtout("",elemstr);
var txtstr = " The Debye temperature allows us to calcuate the specific heat capacity (Eq.6). " 
            + " At room temerature the specific heat capacity is "
            + form(cvrt,1) + " J/kg/K. The graph above shows the temperature dependency of "
            + " the heat capacity. "; 
txtout(txtstr,elemstr);
var txtstr = " The Grüneisen-parameter is set to g=" + form(grueneisen,3) + ". "
            + " We can estimate with the bulk modulus and the heat capacity the linear"
            + " thermal expansoion coefficient with a= " + form(alphal,2) + " µm/m/K."; 
txtout(txtstr,elemstr);



// Fill the table
txtout(compoundname,"c_compoundname");
txtout(form(alat*anull,3)+" Å","c_alat");
txtout(form(volnull,1)+" Å^3","c_volnull");
txtout(strucType,"c_strucType");
txtout(crystForm,"c_crystForm");
txtout(numAtom,"c_numAtom");
txtout(form(dichte,3)+ " kg/l","c_dichte");
txtout(form(molmasse,1)+" g/mol","c_molmasse");
txtout(form(bulkmodul,1)+" GPa","c_bulkmodul");
txtout(form(shearmodul,1)+" GPa","c_shearmodul");
txtout(form(youngmodul,1)+" GPa","c_youngmodul");
txtout(form(c11modul,1)+" GPa","c_c11modul");
txtout(form(hardness,0)+" ","c_hardness");
txtout(form(schall,1)+" m/s","c_schall");
txtout(form(vtwav,1)+" m/S","c_vtwav");
txtout(form(vlwav,1)+" m/S","c_vlwav");
txtout(form(grueneisen,2)+" ","c_grueneisen");
txtout(form(debye,1)+" K","c_debye");
txtout(form(alphal,2)+" µm/m/K","c_alphal");
txtout(form(tschmelz,1)+" K","c_tschmelz");
txtout(form(cvrt,1)+" J/(kg*K)","c_cvrt");
txtout(form(vaperg,1)+" kJ/mol","c_vaperg");
txtout(form(bandwidth,2)+" eV","c_bandwidth");
txtout(form(efermi,2)+" eV","c_efermi");
txtout(form(dosatef,3)+" e¯/eV","c_dosatef");
txtout(form(gapdirect,2)+" eV","c_gapdirect");
txtout(form(nindex,2)+" ","c_nindex");
txtout(form(optigap,2)+" eV","c_optigap");
txtout(form(reflection*100,1)+" %","c_reflection");
txtout(form(transmission*100,1)+" %","c_transmission");
txtout(form(absorpedge,1)+" nm","c_absorpedge");



 //   ffermi=1./(exp((efspace)/(kboltz*tspace(it)))+1)
 //   chgmov(it) =  sum(dosfit.*abs(ffermi-eweight))*defbox + 1e-30 %+ dosatef*defbox % %rhoint*sum(dosfit.*(1-ffermi).*ffermi)*defbox/nvalence 
 //   tauex = 5
 //   phono(it) = ((tspace(it)/debye)^tauex*sum(xdt.^tauex./((exp(xdt)-1).*(1-exp(-xdt))))*dtx)



//eatom = eatom/ehart
//chargeint = sum(qvalnl(:,1)  + qvalnl(:,2).*0.33) %+ qvalnl(:,3).*(10-qvalnl(:,3))/12)
//rhoint = chargeint/volnull
//evalence = sum((efspace+efermi).*dosfit.*eweight)*defbox


//
// plot Bandstructure
var lineChartData = {
    labels: eval(blable),
    datasets: [
        {
            fillColor: "rgba(255, 255, 255, 0)",
            strokeColor: "rgba(000,255,000,1)",
            data: eval(eevk[0])
        }
    ]
}

for (il = 1; il <= lmsum; il++) {
    lineChartData.datasets[il] = {strokeColor: "rgba(000,255,000,1)", "data": eevk[il] };
    if (il > nvalence/2-1) {
        lineChartData.datasets[il] = {strokeColor: "rgba(255,000,000,1)", "data": eevk[il] };
    }
}
var myLine = new Chart(document.getElementById("bandstructure").getContext("2d")).Line(lineChartData);

var barChartData = {
    labels : eval(dlabel),
    datasets : [
        {
            fillColor : "rgba(255,255,255,1)",
            strokeColor : "rgba(45,53,65,1)",
//            highlightFill: "rgba(220,220,220,0.75)",
//            highlightStroke: "rgba(220,220,220,1)",
            data : eval(dosfit)
        }
    ]

}

var myLine = new Chart(document.getElementById("densityofstates").getContext("2d")).Bar(barChartData);



var lineChartData = {
    labels: eval(tlabel),
    datasets: [
        {
            fillColor: "rgba(255, 255, 255, 0)",
            strokeColor: "rgba(000,255000,1)",
            data: eval(cvt)
        }
    ]
}
var myLine = new Chart(document.getElementById("debyetemp").getContext("2d")).Line(lineChartData);

var lineChartData = {
    labels: eval(olabel),
    datasets: [
        {
            fillColor: "rgba(255, 255, 255, 0)",
            strokeColor: "rgba(000,000,255,1)",
            data: eval(pjdos)
        }
    ]
}
var myLine = new Chart(document.getElementById("pjdosplot").getContext("2d")).Line(lineChartData);

var lineChartData = {
    labels: eval(lambda),
    datasets: [
        {
            fillColor: "rgba(255, 255, 255, 0)",
            strokeColor: "rgba(255,000,255,1)",
            data: eval(normeps1)
        },
        {
            fillColor: "rgba(255, 255, 255, 0)",
            strokeColor: "rgba(000,255,255,1)",
            data: eval(normeps2)
        }
    ]
}
var myLine = new Chart(document.getElementById("eps1eps2").getContext("2d")).Line(lineChartData);

// References
// https://www.researchgate.net/publication/280131674_Simple_Tight_Binding_Electronic_Structure_Calculation_Method_based_on_a_Mathematical_Prototyping_System
// https://en.wikipedia.org/wiki/Debye_model
// https://en.wikipedia.org/wiki/Tight_binding
// https://en.wikipedia.org/wiki/Kramers%E2%80%93Kronig_relations
// https://en.wikipedia.org/wiki/Refractive_index  siehe Relative permittivity and permeability
// https://en.wikipedia.org/wiki/Bulk_modulus
// https://en.wikipedia.org/wiki/Anton-Schmidt_equation_of_state

}