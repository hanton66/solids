//
// Some genral Math funktions for easier typing
//
// Square-Root
function sqrt(x) {
    var y;
    y = Math.sqrt(x);
    return y;
}

//
// Power-Function
function power(x, z) {
    var y;
    y = Math.pow(x, z);
    return y;
}

//
// Exp-Funktion
function exp(x) {
    var y;
    y = Math.exp(x);
    return y;
}
//
// log-Funktion
function log(x) {
    var y;
    y = Math.log(x);
    return y;
}
//
// complex exponteial funktion y=exp(i*x)
function expi(x) {
var y=[];
y[1] = Math.cos(x);
y[2] = Math.sin(x)
return y;
}

//
// Sinus-Hyberbolicus-Function
function sinh(x) {
    var y;
    y = 0.5 * (exp(x) - exp(-x));
    return y;
}
//
// Tangens-Hyberbolicus-Function
function tanh(x) {
    var y;
    y = (exp(2*x) - 1)/(exp(2*x) + 1);
    return y;
}

//
// Nummeric Differential
function diff(fofr, r, ngo, cmulti) {
    var dfdr = new Array(ngo);
    var dr = 0;
    for (igo = 1; igo < ngo; igo++) {
        dr = r[igo+1] - r[igo-1];
        dfdr[igo] = (fofr[igo+1] - fofr[igo-1])/dr*cmulti;
    }
    dfdr[0] = 2*dfdr[1]-dfdr[2];
    dfdr[ngo] = 2*dfdr[ngo-1]-dfdr[ngo-2];
    return dfdr;
}

//
// Nummeric Second Derivation
function diff2(fofr, r, ngo, cmulti) {
    var d2fdr2 = new Array(ngo);
    var dfdr = new Array(ngo);    
    dfdr = diff(fofr, r, ngo, cmulti);
    d2fdr2 = diff(dfdr, r, ngo, cmulti);
    return d2fdr2;
}

function chkCross(fofr,r,refval) {
    var rcross = 0.0;
    var dr = 0.0;
    var dall = 0.0;
    var dfofr = 0.0;
    var ngo = fofr.length;
    for (igo = 0; igo < ngo; igo++) {
        if (fofr[igo] < refval) {
            if (fofr[igo+1] > refval) {
                var dr = r[igo+1] - r[igo]; 
                var dall = fofr[igo+1] - fofr[igo];
                var dfofr = refval - fofr[igo]; 
                rcross = r[igo] + dr*dfofr/dall;
            }
        }
        if (fofr[igo] > refval) {
            if (fofr[igo+1] < refval) {
                var dr = r[igo+1] - r[igo]; 
                var dall = fofr[igo+1] - fofr[igo];
                var dfofr = refval - fofr[igo]; 
                rcross = r[igo] + dr*dfofr/dall;
            }
        }     
        if (fofr[igo] == refval) {
            rcross = r[igo];
        }   
    }
    return rcross;

}

//
//  chk same value
function isequal(x,y) {
    var z = 0;
    if (x == y) {
        z = 1;
    }
    return z;
}

//
// Sort an one-dimensional array
function sortV(x_in) {
    var y = []; x = x_in;
    var lt=x.length;
    for (jgo = 1; jgo < lt; jgo++) {
        for (igo = 1; igo < lt; igo ++) {
            if (x[igo] < x[igo-1]) {
                tmp = x[igo];
                x[igo] = x[igo-1];
                x[igo-1] = tmp;
            }
        }
    }
    y = x;
    return y;
    }
//
// get max Valueof one-dimensional array
function maxV(x) {
    var xerg = x[0]; y = xerg;
    var lt=x.length;
    for (igo = 1; igo < lt; igo++) {
        if (x[igo] > xerg) {
            xerg = x[igo]; 
        }
    }
    y = xerg;
    return y;
}//
// get min Valueof one-dimensional array    
function minV(x) {
    var xerg = x[0]; y = xerg;
    var lt=x.length;
    for (igo = 1; igo < lt; igo++) {
        if (x[igo] < xerg) {
            xerg = x[igo]; 
        }
    }
    y = xerg;
    return y;
}  
//
// get the sum of all values of vector
function sumV(x) {
    var xsum = 0;
    var lt=x.length;
    for (igo = 0; igo < lt; igo++) {
        xsum = xsum + x[igo];
    }
    return xsum;
}  


//
// Create a histogramm
function hist(x,n) {
    var y = numeric.rep([n],0);
    var lt=x.length;
    var xmax = Math.ceil(maxV(x)/5)*5;  
    var xmin = Math.floor(minV(x)/5)*5; 
    var dbox = (xmax-xmin)/n;
    var ibox = 1;
    for (igo = 0; igo < lt; igo++) {
        ibox = Math.floor((x[igo]-xmin)/dbox);
        y[ibox] = y[ibox] + 1;
    }
    return y;
}



//
// get the max values of an array
function maxM(x_in) {
    var y = 0; 
    var lt=x_in.length;
    var xmin = x_in[0][0];
    for (igo = 0; igo < lt; igo++) {
        for (jgo = 0; jgo < lt; jgo ++) {
            if (x_in[igo][jgo] > xmin) {
               xmin = x_in[igo][jgo];
            }
        }
    }
    y = xmin;
    return y;
}

//
// get the sum of all absolute values of matrix
function sumaM(x_in) {
    var y = 0; 
    var lt=x_in.length;
    var y = 0.0;
    for (igo = 0; igo < lt; igo++) {
        for (jgo = 0; jgo < lt; jgo ++) {
            y = y + Math.abs(x_in[igo][jgo])
        }
    }
    return y;
}

//
// Check symmetrie of the matrix
function chkMatSy(x_in) {
    var y = 0; 
    var lt=x_in.length;
    var y = 0.0;
    for (igo = 0; igo < lt; igo++) {
        for (jgo = igo+1; jgo < lt; jgo ++) {
            y = y + Math.abs(x_in[igo][jgo] + x_in[jgo][igo]); 
        }
    }
    return y;
}


//
// Formatting Number
function form(x, z) {
    var y;
    x = Math.round(x * power(10, z));
    y = x / power(10, z);
    return y;
}

//
//
function txtout(txtstr,elemstr) {
var nxtline = document.createElement('br');
var element = document.getElementById(elemstr).appendChild(nxtline); 
var wrtline = document.createTextNode(txtstr);
var element = document.getElementById(elemstr).appendChild(wrtline);
}