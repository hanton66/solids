function addAtom(getAtom) {

    noAtom = noAtom + 1;

    if (noAtom <= maxCompund) {

    var nxtline = document.createElement('br');
    var wrtline = document.createTextNode("Added " + getAtom + " to the structure.");
    var element = document.getElementById('CrystalAtoms').appendChild(nxtline); 
    var element = document.getElementById('CrystalAtoms').appendChild(wrtline);  
    
    if (getAtom == "H") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 1,
                            "e_spd": [-0.233 ,  0.000,  0.000],
                            "q_spd": [     1,       0,      0], 
                            "ratom":  31, 
                            "nexp": 2.0 };  
    }
 
    if (getAtom == "He") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 2,
                            "e_spd": [-0.570 ,  0.000,  0.000],
                            "q_spd": [     2,       0,      0], 
                            "ratom":  50, 
                            "nexp": 2.0 };   
    }  

    if (getAtom == "Li") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 3,
                            "e_spd": [-0.106 ,  0.000,  0.000],
                            "q_spd": [     1,       0,      0], 
                            "ratom": 152, 
                            "nexp": 1.7 };  
    }

    if (getAtom == "Be") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 4,
                            "e_spd": [-0.206 ,  0.000,  0.000],
                            "q_spd": [     2,       0,      0], 
                            "ratom": 111, 
                            "nexp": 1.6 };  
    }    
  
    if (getAtom == "B") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 5,
                            "e_spd": [-0.344, -0.137,  0.000],
                            "q_spd": [     2,       1,      0], 
                            "ratom":  84, 
                            "nexp": 1.6 };  
    }
    

    if (getAtom == "C") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 6,
                            "e_spd": [-0.500, -0.199,  0.000],
                            "q_spd": [     2,       2,      0], 
                            "ratom":  76, 
                            "nexp": 1.7 };  
    }


    if (getAtom == "N") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 7,
                            "e_spd": [-0.676, -0.266,  0.000],
                            "q_spd": [     2,       3,      0], 
                            "ratom":  71, 
                            "nexp": 2.0 };  
    }
    
    if (getAtom == "O") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 8,
                            "e_spd": [-0.871 , -0.338,  0.000],
                            "q_spd": [     2,       4,      0], 
                            "ratom":  66, 
                            "nexp": 2.0 };  
    }     

    if (getAtom == "F") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 9,
                            "e_spd": [-1.087, -0.416,  0.000],
                            "q_spd": [     2,       5,      0], 
                            "ratom":  57, 
                            "nexp": 1.5 };  
    }

    if (getAtom == "Ne") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 3,
                            "e_spd": [-1.323 , -0.498,  0.000],
                            "q_spd": [     2,       6,      0], 
                            "ratom":  50, 
                            "nexp": 2.0 };  
    }

    if (getAtom == "Na") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 11,
                            "e_spd": [-0.103, -0.000,  0.000],
                            "q_spd": [     1,       0,      0], 
                            "ratom":  186, 
                            "nexp": 1.9 }; 
                        }
  
    if (getAtom == "Ca") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 12,
                            "e_spd": [-0.175, -0.000,  0.000],
                            "q_spd": [     2,       0,      0], 
                            "ratom":  160, 
                            "nexp": 2.0 }; 
                        }


    if (getAtom == "Al") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 13,
                            "e_spd": [-0.287, -0.103,  0.000],
                            "q_spd": [     2,       1,      0], 
                            "ratom":  143, 
                            "nexp": 1.85 };     
                        }


    if (getAtom == "Si") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 14,
                            "e_spd": [-0.398, -0.153,  0.000],
                            "q_spd": [     2,       2,      0], 
                            "ratom":  111, 
                            "nexp": 1.9 };  
    }

    if (getAtom == "P") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 15,
                            "e_spd": [-0.512, -0.206,  0.000],
                            "q_spd": [     2,       3,      0], 
                            "ratom":  107, 
                            "nexp": 2.0 };      
                        }


    if (getAtom == "S") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 16,
                            "e_spd": [-0.631, -0.261,  0.000],
                            "q_spd": [     2,       4,      0], 
                            "ratom":  105, 
                            "nexp": 2.0 };      
                        }

    if (getAtom == "Cl") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 17,
                            "e_spd": [-0.754, -0.320,  0.000],
                            "q_spd": [     2,       5,      0], 
                            "ratom":  75, 
                            "nexp": 1.4 };     
                        }


    if (getAtom == "Ga") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 31,
                            "e_spd": [-0.328, -0.102,  0.000],
                            "q_spd": [     2,       1,      0], 
                            "ratom":  122, 
                            "nexp": 2.0 };  
    }
    

 
    if (getAtom == "Ge") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 32,
                            "e_spd": [-0.427, -0.150,  0.000],
                            "q_spd": [     2,       2,      0], 
                            "ratom":  122, 
                            "nexp": 2.1 };  
    } 

    if (getAtom == "As") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 32,
                            "e_spd": [-0.523, -0.197,  0.000],
                            "q_spd": [     2,       3,      0], 
                            "ratom":  119, 
                            "nexp": 2.0 };  
    }


    }

    if (noAtom > maxCompund) {
        noAtom = noAtom - 1;
        var nxtline = document.createElement('br');
        var wrtline = document.createTextNode("Maximal number of Atoms reached.");
        var element = document.getElementById('CrystalAtoms').appendChild(nxtline); 
        var element = document.getElementById('CrystalAtoms').appendChild(wrtline);  
    }

    if (noAtom == maxCompund) {
        document.getElementById('readyReport').disabled = false; 
    }

    console.log (espdofz);
}

function delAtom() {
    if (noAtom > 0) {
    atomSort = espdofz[noAtom].atName
    noAtom = noAtom - 1;
    var nxtline = document.createElement('br');
    var wrtline = document.createTextNode("Deleted the last Atom " + atomSort + ".");
    var element = document.getElementById('CrystalAtoms').appendChild(nxtline); 
    var element = document.getElementById('CrystalAtoms').appendChild(wrtline);
    }
}