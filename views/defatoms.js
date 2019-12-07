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
                            "ratom":  250, 
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
                            "ratom":  220, 
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
  
    if (getAtom == "Mg") {
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


    if (getAtom == "Ar") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 18,
                            "e_spd": [-0.883, -0.382,  0.000],
                            "q_spd": [     2,       6,      0], 
                            "ratom":  106, 
                            "nexp": 2.0 };     
                        }

    
    if (getAtom == "K") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 19,
                            "e_spd": [-0.089, - 0.000,  0.000],
                            "q_spd": [     1,       0,      0], 
                            "ratom":  227, 
                            "nexp": 2.1 };     
                           }


    if (getAtom == "Ca") {
          espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 20,
                            "e_spd": [-0.141,  0.000,  0.000],
                            "q_spd": [     2,       0,      0], 
                            "ratom":  197, 
                            "nexp": 2.1 };     
                            }

 
    if (getAtom == "Sc") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 21,
                            "e_spd": [-0.156,  0.000,  -0.131],
                            "q_spd": [     2,       0,      1], 
                            "ratom":  161, 
                            "nexp": 2.3 };     
                         }

    if (getAtom == "Ti") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 22,
                            "e_spd": [-0.167, -0.000,  -0.170],
                            "q_spd": [     2,       0,      2], 
                            "ratom":  145, 
                            "nexp": 2.2 };     
                        }                        
                        


    if (getAtom == "V") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 23,
                            "e_spd": [-0.176,  0.000,  -0.204],
                            "q_spd": [     2,       0,      3], 
                            "ratom":  131, 
                            "nexp": 2.2 };     
                        }


    if (getAtom == "Cr") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 24,
                            "e_spd": [-0.150, -0.000,  -0.118],
                            "q_spd": [     2,       0,      4], 
                            "ratom":  124, 
                            "nexp": 2.3 };     
                        }

    if (getAtom == "Mn") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 25,
                            "e_spd": [-0.191, -0.000,  -0.267],
                            "q_spd": [     2,       0,      5], 
                            "ratom":  137, 
                            "nexp": 2.2 };     
                        }

 
    if (getAtom == "Fe") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 26,
                            "e_spd": [-0.198,  0.000,  -0.296],
                            "q_spd": [     2,       0,      6], 
                            "ratom":  124, 
                            "nexp": 2.1 };     
                        }


    if (getAtom == "Co") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 27,
                            "e_spd": [-0.204,  0.000,  -0.322],
                            "q_spd": [     2,       0,      7], 
                            "ratom":  125, 
                            "nexp": 2.0 };     
                        }

  
    if (getAtom == "Ni") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 28,
                            "e_spd": [-0.210,  0.000,  -0.349],
                            "q_spd": [     2,       0,      8], 
                            "ratom":  125, 
                            "nexp": 2.0 };     
                        }


    if (getAtom == "Cu") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 29,
                            "e_spd": [-0.172, 0.000,  -0.202],
                        "q_spd": [     1,       0,     10], 
                        "ratom":  128, 
                                        "nexp": 2.0 };  
                        }  


 
    if (getAtom == "Zn") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 30,
                            "e_spd": [-0.222,  0.000,  -0.399],
                            "q_spd": [     2,       0,     10], 
                            "ratom":  133, 
                            "nexp": 1.8 };     
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
                            "atCharge": 33,
                            "e_spd": [-0.523, -0.197,  0.000],
                            "q_spd": [     2,       3,      0], 
                            "ratom":  119, 
                            "nexp": 2.0 };  
    }

    if (getAtom == "Au") {
        espdofz[noAtom] = { "atName": getAtom,
                            "atCharge": 79,
                            "e_spd": [-0.162,  0.000,  -0.306],
                            "q_spd": [     1,       0,     10], 
                            "ratom":  144, 
                            "nexp": 2.5 };  
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