/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.

*/

/**
 * JSXGraph namespace. Holds all classes, objects, functions and variables belonging to JSXGraph
 * to reduce the risc of interfering with other JavaScript code.
 * @namespace
 */
var JXG = {};
(function(){
    var i, s;
    //JXG.useMinify = true;
    JXG.countDrawings = 0;
    JXG.countTime = 0;
    JXG.require = function(libraryName) {};
    JXG.rendererFiles = [];
    JXG.rendererFiles['svg'] = 'SVGRenderer';
    JXG.rendererFiles['vml'] = 'VMLRenderer';
    JXG.baseFiles = null;
    // this maybe required by additional software/extensions and/or future renderers
    JXG.requirePath = '';
    for (i=0;i<document.getElementsByTagName("script").length;i++) {
        s = document.getElementsByTagName("script")[i];
        if (s.src && s.src.match(/loadjsxgraphInOneFile\.js(\?.*)?$/)) {
            JXG.requirePath = s.src.replace(/loadjsxgraphInOneFile\.js(\?.*)?$/,'');
        }
    }
JXG.serverBase = JXG.requirePath + 'server/';
})();
/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph. If not, see <http://www.gnu.org/licenses/>.
*/

/** 
 * @fileoverview In this file the namespace JXG.Math is defined, which is the base namespace
 * for namespaces like Math.Numerics, Math.Algebra, Math.Statistics etc.
 * @author graphjs
 */
 
 /**
  * Math namespace.
  */
JXG.Math = new Object();

/* Math constants */
JXG.Math.eps = 0.000001;

/**
 * Represents a vector.
 * @constructor
 * @param {Array} elements An array of numerical values containing the coefficients to be put into the vector.
 * @return JXG.Math.Vector
 */
JXG.Math.Vector = function(elements) {
    var i;
    this.length = 0;
    
    if((typeof elements != undefined) && (elements != null)) {
        for(i=0; i<elements.length; i++) {
            this.push(elements[i]);
        }
    }
};

/*
 * The base class for Vector is just an array.
 */
JXG.Math.Vector.prototype = new Array();

/**
 * Returns the dimension of the vector.
 * @type int
 */
JXG.Math.Vector.prototype.n = function() {
    return this.length;
};

/**
 * Exchanges two elements of the vector.
 * @param {int} i The first element that is to be exchanged.
 * @param {int} j The second element that is to be exchanged.
 */
JXG.Math.Vector.prototype.exchange = function(i, j) {
    var temp = this[i];
    
    this[i] = this[j];
    this[j] = temp; 
};

/**
 * Represents a matrix.
 * @constructor
 * @param {Array} elements An 2-dimensional array of numerical values containing the coefficients to be put into the vector.
 * @throws {JXG.DimensionMismatchException} If the rows of the matrix don't have all the same length.
 * @return JXG.Math.Vector
 */
JXG.Math.Matrix = function(elements) {
    var oldLength = 0,
        testLength = false,
        i, j, len, leni;
    
    this.length = 0;
    
    if ((typeof elements != undefined) && (elements != null)) {
        len = elements.length;
        for (i=0; i<len; i++) {
            leni = elements[i].length;
            this.push(new Array());
          
            if (testLength) {
                if (oldLength != leni) {
                    this.length = 0;
                    throw new JXG.DimensionMismatchException("Your array contains arrays with different lengths.");
                }
            }
                
            for (j=0; j<leni; j++) {
                this[i].push(elements[i][j]);
            }
          
            oldLength = leni;
            testLength = true;
        }
    }
};

/*
 * The base class for Matrix is also just an array.
 */
JXG.Math.Matrix.prototype = new Array();

/**
 * Returns the amount of rows of the matrix.
 * @type int
 */
JXG.Math.Matrix.prototype.m = function() {
    return this.length;
};

/**
 * Returns the amount of columns of the matrix.
 * @type int
 */
JXG.Math.Matrix.prototype.n = function() {
    if(this.length > 0)
        return this[0].length;
    else
        return 0;
};

/**
 * Exchanges two rows of the matrix.
 * @param {int} i The first row that is to be exchanged.
 * @param {int} j The second row that is to be exchanged.
 */
JXG.Math.Matrix.prototype.exchangeRows = function(i, j) {
   var temp = this[i];
    
   this[i] = this[j];
   this[j] = temp; 
};

/**
 * Exception signaling inconsistent dimension conditions.
 * @constructor
 * @param {string} message A message which explains what went wrong-
 */
JXG.DimensionMismatchException = function(message) {
    if ((typeof message != undefined) && (message != null))
        this.message = message;
    else
        this.message = null;
};

/**
 * Returns a string explaining, what exactly went wrong.
 * @type string
 * @return A string explaining why this exception was raised.
 */
JXG.DimensionMismatchException.prototype.what = function() {
    var default_msg = "Matrix has incorrect dimensions";
   
    if (this.message != null)
        return default_msg + ": " + this.message + ".";
    else
        return default_msg + ".";
};

/**
 * Exception signaling an singular matrix.
 * @constructor
 * @param {string} message A message which explains what exactly went wrong-
 */
JXG.SingularMatrixException = function(message) {
    if ((typeof message != undefined) && (message != null))
        this.message = message;
    else
        this.message = null;
};

/**
 * Returns a string explaining, what exactly went wrong.
 * @type string
 * @return A string explaining why this exception was raised.
 */
JXG.SingularMatrixException.prototype.what = function() {
    var default_msg = "Matrix is singular";
   
    if (this.message != null)
        return default_msg + ": " + this.message + ".";
    else
        return default_msg + ".";
};


/**
 * Matrix-vector multiplication.
 * @param {Array} mat1 Two dimensional array of numbers
 * @param {Array} vec Array of numbers
 * @return {Array} Array of numbers containing result
 */
JXG.Math.matVecMult = function(/** array */ mat1, /** array */ vec) /** array */ {
    var m = mat1.length,
        n = vec.length,
        res = [],
        i, s, k;
    if (n==3) {
        for (i=0;i<m;i++) {
            res[i] = mat1[i][0]*vec[0] + mat1[i][1]*vec[1] + mat1[i][2]*vec[2];
        }
    } else {
        for (i=0;i<m;i++) {
            s = 0;
            for (k=0;k<n;k++) { s += mat1[i][k]*vec[k]; }
            res[i] = s;
        }
    }
    return res;
};

/**
 * Matrix-matrix multiplication.
 * @param {Array} mat1 Two dimensional array of numbers
 * @param {Array} mat2 Two dimensional array of numbers
 * @return {Array} Two dimensional Array of numbers containing result
 */
JXG.Math.matMatMult = function(/** array */ mat1, /** array */ mat2) /** array */ {
    var m = mat1.length,
        n = mat2[0].length,
        m2 = mat2.length,
        res = [], 
        i, j, s, k;
        
    for (i=0;i<m;i++) {
        res[i] = [];
    }

    for (i=0;i<m;i++) {
        for (j=0;j<n;j++) {
            s = 0;
            for (k=0;k<m2;k++) {
                s += mat1[i][k]*mat2[k][j];
            }
            res[i][j] = s;
        }
    }
    return res;
};

/**
 * Transpose a matrix which is of type array of arrays.
 * @param {Array} M 
 * @return {Array} transpose of M
 */
JXG.Math.Matrix.transpose = function(/** Array */ M) /** Array*/  {
    var MT = [], i, j, 
        m, n;
    
    m = M.length;                   // number of rows of M
    n = (M.length>0)?M[0].length:0; // number of columns of M

    for (i=0;i<n;i++) {
        MT.push([]);
        for (j=0;j<m;j++) {
            MT[i].push(M[j][i]);
        }
    }
    return MT;
}

/**
  * Calculates the crossproducts of two vectors
  * of length three.
  * In case of homogeneous coordinates this is either
  * - the intersection of two lines
  * - the line through two points.
  * @param {Array} c1 homogeneous coordinates of line (point) 1
  * @param {Array} c2 homogeneous coordinates of line (point) 2
  * @type Array
  * @return vector of length 3:  homogeneous coordinates
  *   of the resulting line / point.
  */
JXG.Math.crossProduct = function(c1,c2) {
    return [c1[1]*c2[2]-c1[2]*c2[1],
            c1[2]*c2[0]-c1[0]*c2[2],
            c1[0]*c2[1]-c1[1]*c2[0]];
};

/**
 * Inner product of two vectors a, b. n is the length of the vectors.
 * @param a Vector
 * @param b Vector
 * @param [n] Length of the Vectors. If not given the length of the first vector is taken.
 * @return The inner product of a and b.
 */
JXG.Math.innerProduct = function(a, b, n) {    
    var i, s = 0;
    
    if(typeof n == 'undefined')
        n = a.length;
    
    for (i=0;i<n;i++) {
        s += a[i]*b[i];
    }
    return s;
};



/**
* Dynamic programming approach for recursive functions.
* From "Speed up your JavaScript, Part 3" by Nicholas C. Zakas.
* @see JXG.Math.factorial
* http://blog.thejit.org/2008/09/05/memoization-in-javascript/
*/
JXG.memoizer = function (f) {
    var cache, join;
    
    if (f.memo) {
        return f.memo;
    }
    cache = {};
    join = Array.prototype.join;

    return (f.memo = function() {
        var key = join.call(arguments);
        //return (key in cache)
        return (typeof cache[key]!='undefined') // Seems to be a bit faster than "if (a in b)"
            ? cache[key]
            : cache[key] = f.apply(this, arguments);
    });
};

/**
* Compute the factorial of a positive integer.
* @param {integer n}
* @return {return n*(n-1)...2*1}
*/
JXG.Math.factorial = JXG.memoizer(function (n) {
        if (n<0) return NaN; 
        if (n==0 || n==1) return 1;
        return n*arguments.callee(n-1);
});

/**
* Comupte the binomial coefficient.
* @param {integer n}
* @param {integer k}
* 
* @return {n\choose k}
*/
JXG.Math.binomial = JXG.memoizer(function(n,k) {
    var b, i;
    
    if (k>n || k<0) return 0;
    if (k==0 || k==n) return 1;
    
    b = 1;
    for (i=0;i<k;i++) {
        b *= (n-i);
        b /= (i+1);
    }
    return b;
    //return arguments.callee(n-1,k-1)+arguments.callee(n-1,k);
});

/*
    // Just for test purposes;
    
JXG.Math.Numerics.prototype.fibonacci = JXG.memoizer(function (n) {
        if(n < 2) return 1; else return arguments.callee(n-2) + arguments.callee(n-1);  
    });
*/    

/**
* Round a decimal number to n decimal places
* @deprecated Use (number).toFixed(n) instead.
* @param {float num} Number to round
* @param {integer n} number of digits after the point to leave
* 
* @return {rounded num}
*/
JXG.Math.round = function(num, n) {
    var z, s;
    //return Math.round(num*Math.pow(10,n))/Math.pow(10,n);
    //var z = num.toFixed(n);
    
    z = num - Math.ceil(num);
    s = z.toString();
    if (z < 0) {
        s = s.substr(0,n+3);
    }
    else {
        s = s.substr(0,n+2);
    }
    z = parseFloat(s);
    t = parseInt(num.toString());
    return t+z;
};

/**
 * Cosine hyperbolicus of x.
 * @param {float} x The number the cosine hyperbolicus will be calculated of.
 * @return {float} Cosine hyperbolicus of the given value.
 */
JXG.Math.cosh = function(/** number */ x) /** number */ {
    return (Math.exp(x)+Math.exp(-x))*0.5;
};

/**
 * Sine hyperbolicus of x.
 * @param {number} x The number the sine hyperbolicus will be calculated of.
 * @return {number} Sine hyperbolicus of the given value.
 */
JXG.Math.sinh = function(/** number */ x) /** number */ {
    return (Math.exp(x)-Math.exp(-x))*0.5;
};


/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph. If not, see <http://www.gnu.org/licenses/>.
*/

/** 
 * @fileoverview In this file the namespace Math.Numerics is defined, which holds numerical
 * algorithms for solving linear equations etc.
 * @author graphjs
 */


/**
 * Math.Numerics namespace holds numerical algorithms, constants, and variables.
 * @namespace
 */
JXG.Math.Numerics = {}; 

/**
 * A constant representing the trapez rule for numerical integration. Is used in conjunction with {@link JXG.Math.Numerics.integration_type}
 * and {@link JXG.Math.Numerics.NewtonCotes}.
 * @type number
 * @constant
 */
JXG.Math.Numerics.INT_TRAPEZ  = 0x00001;

/**
 * A constant representing the simpson rule for numerical integration. Is used in conjunction with {@link JXG.Math.Numerics.integration_type}
 * and {@link JXG.Math.Numerics.NewtonCotes}.
 * @type number
 * @constant
 */
JXG.Math.Numerics.INT_SIMPSON = 0x00002;

/**
 * A constant representing the milne rule for numerical integration. Is used in conjunction with {@link JXG.Math.Numerics.integration_type}
 * and {@link JXG.Math.Numerics.NewtonCotes}.
 * @type number
 * @constant
 */
JXG.Math.Numerics.INT_MILNE   = 0x00003;
  
/**
 * Number of nodes for evaluation, used for integration algorithms.
 * @type number
 */
JXG.Math.Numerics.number_of_nodes = 28;

/**
 * Type of integration algorithm, possible values are:
 * <ul>
 *   <li>{@link JXG.Math.Numerics.INT_TRAPEZ}</li>
 *   <li>{@link JXG.Math.Numerics.INT_SIMPSON}</li>
 *   <li>{@link JXG.Math.Numerics.INT_MILNE}</li>
 * </ul>
 * @type number
 */
JXG.Math.Numerics.integration_type = JXG.INT_MILNE;

/**
 * Solves a system of linear equations given by the right triangular matrix R and vector b.
 * @param R Right triangular matrix. All entries a_(i,j) with i < j are ignored.
 * @param b Right hand side of the linear equation system.
 * @return A vector that solves the system of linear equations.
 * @private
 */ 
JXG.Math.Numerics.backwardSolve = function(/** JXG.Math.Matrix */ R, /** JXG.Math.Vector */ b) /** JXG.Math.Vector */ {
    var x = b,
        m, n,
        i, j;

    // m: number of rows of R
    // n: number of columns of R
    // Relaxation: R may be of type JXG.Math.Matrix or Array
    //             b may be of type JXG.Math.Vector or Array
    if (R.m) { // R is of type JXG.Math.Matrix
        m = R.m();
        n = R.n();
    } else {   // R is of type array
        m = R.length;
        n = (R.length>0)?R[0].length:0;
    }
    for (i = m-1; i >= 0; i--) {
        for (j = n-1; j > i; j--) {
            x[i] -= R[i][j] * x[j];
        }
        x[i] /= R[i][i];
    }
   
    return x;
};

/**
 * Solves a system of linear equations given by A and b using the Gauss-Jordan-elimination.
 * The algorithm runs in-place. I.e. the entries of A and b are changed.
 * @param A Square matrix containing the coefficients of the lineare equation system.
 * @param b A vector containing the linear equation system's right hand side. 
 * @throws {JXG.DimensionMismatchException} If a non-square-matrix is given or the b has not the right length.
 * @throws {JXG.SingularMatrixException} If A's rank is not full.
 * @return A vector that solves the linear equation system.
 */
JXG.Math.Numerics.Gauss = function(/** JXG.Math.Matrix */ A, /** JXG.Math.Vector */ b) /** JXG.Math.Vector */ {
    var eps = JXG.Math.eps,
        n,                 
        i, j, k, P,
        x, y;
    
    // n: number of columns of A
    // Relaxation: A may be of type JXG.Math.Matrix or Array
    //             b may be of type JXG.Math.Vector or Array
    if (A.n) { // A is of type JXG.Math.Matrix
        n = A.n();
    } else {   // A is of type array
        n = (A.length>0)?A[0].length:0;
    }
    
    /* vector to keep track of permutations caused by pivotion */
    P = new JXG.Math.Vector();
    for (i = 0; i < n; i++) {
        P.push(i);
    }
    
   /* Gauss-Jordan-elimination */
    for (j=0; j < n; j++)
    {
        for (i = n-1; i > j; i--) {
            /* Is the element which is to eliminate greater than zero? */
            if (Math.abs(A[i][j]) > JXG.Math.eps) {
                /* Equals pivot element zero? */
                if (Math.abs(A[j][j]) < JXG.Math.eps) {
                    /* Yeah, so we have to exchange the rows */
                    A.exchangeRows(i, j);
                    b.exchange(i, j);
                    P.exchange(i, j);
                }
                else {
                    /* Saves the L matrix of the LR-decomposition. unneeded. */
                    A[i][j] /= A[j][j];
                    /* Transform right-hand-side b */
                    b[i] -= A[i][j] * b[j];
                    /* subtract the multiple of A[i][j] / A[j][j] of the j-th row from the i-th. */
                    for (k = j + 1; k < n; k ++) {
                        A[i][k] -= A[i][j] * A[j][k];
                    }
                }
            }
            if (Math.abs(A[j][j]) < JXG.Math.eps) { // The absolute values of all coefficients below the j-th row in the j-th column are smaller than JXG.Math.eps.
                throw new SingularMatrixException();
            }
        }
    }
   
    return JXG.Math.Numerics.backwardSolve(A, b); // return Array
    /* 
    y = JXG.Math.Numerics.backwardSolve(A, b);
    x = new JXG.Math.Vector();
    for (i = 0; i < n ; i++) { // y.n()
        x.push(y[P[i]]);
    }
   
    return x; // return JXG.Math.Vector
    */
};

/**
 * Compute the inverse of an nxn matrix with Gauss elimination.
 * @param {JXG.Math.Matrix} Ain
 */
JXG.Math.Numerics.Inverse = function(Ain) {
    var i,j,k,s,ma,r,swp,
        n = Ain.length,
        A = [],
        p = [],
        hv = [];
        
    for (i=0;i<n;i++) {
        A[i] = [];
        for (j=0;j<n;j++) { A[i][j] = Ain[i][j]; }
        p[i] = i;
    }
    
    for (j=0;j<n;j++) {
        // pivot search:
        ma = Math.abs(A[j][j]);
        r = j;
        for (i=j+1;i<n;i++) {
            if (Math.abs(A[i][j])>ma) {
                ma = Math.abs(A[i][j]);
                r = i;
            }
        }
        if (ma<=JXG.Math.eps) { // Singular matrix
            return false;  
        }
        // swap rows:
        if (r>j) {
            for (k=0;k<n;k++) { 
                swp = A[j][k]; A[j][k] = A[r][k]; A[r][k] = swp;
            }
            swp = p[j]; p[j] = p[r]; p[r] = swp;
        }
        // transformation:
        s = 1.0/A[j][j];
        for (i=0;i<n;i++) {
            A[i][j] *= s;
        }
        A[j][j] = s;
        for (k=0;k<n;k++) if (k!=j) {
            for (i=0;i<n;i++) if (i!=j) {
                A[i][k] -= A[i][j]*A[j][k];
            }
            A[j][k] = -s*A[j][k];
        }
    }
    // swap columns:
    for (i=0;i<n;i++) {
        for (k=0;k<n;k++) { hv[p[k]] = A[i][k]; }
        for (k=0;k<n;k++) { A[i][k] = hv[k]; }
    }
    return A;
};

/**
 * NEEDS IMPLEMENTATION. TODO Decomposites the matrix A in an orthogonal matrix Q and a right triangular matrix R. 
 * @param {JXG.Math.Matrix} A A matrix.
 * @type Object
 * @throws {Exception} If A's rank is not full.
 * @return The matrices Q and R.
 * @private
 */
JXG.Math.Numerics.QR = function(A, b) {
    // TODO needs implementation
};

/**
 * Compute the Eigenvalues and Eigenvectors of a symmetric 3x3 matrix with the Jacobi method
 * @param {JXG.Math.Matrix} Ain A symmetric 3x3 matrix.
 * @type Object
 * @throws {Exception} If A's rank is not full.
 * @return [A,V] the matrices A and V. The diagonal of A contains the Eigenvalues,
 * V contains the Eigenvectors.
 */
JXG.Math.Numerics.Jacobi = function(Ain) {
    var i,j,k,ih,aa,si,co,tt,
        sum = 0.0,
        ssum, amax,
        n = Ain.length,
        V = [[0,0,0],[0,0,0],[0,0,0]],
        A = [[0,0,0],[0,0,0],[0,0,0]];
    
    // Initialization. Set initial Eigenvectors.
    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) { 
            V[i][j] = 0.0;
            A[i][j] = Ain[i][j];
            sum += Math.abs(A[i][j]);
        }
        V[i][i] = 1.0;
    }
    // Trivial problems
    if (n==1) { return [A,V]; }
    if (sum<=0.0) { return [A,V]; }
    
    sum /= (n*n);

    // Reduce matrix to diagonal
    do {
        ssum = 0.0;
        amax = 0.0;
        for (j=1;j<n;j++) {
            for (i=0;i<j;i++) {
                // Check if A[i][j] is to be reduced
                aa = Math.abs(A[i][j]);
                if (aa>amax) { amax = aa; }
                ssum += aa;
                if (aa<0.1*amax) {
                    continue;
                } else {
                    // calculate rotation angle
                    aa = Math.atan2(2.0*A[i][j],A[i][i]-A[j][j])*0.5;
                    si = Math.sin(aa);
                    co = Math.cos(aa);
                    // Modify 'i' and 'j' columns
                    for (k=0;k<n;k++) {
                        tt = A[k][i];
                        A[k][i] =  co*tt+si*A[k][j];
                        A[k][j] = -si*tt+co*A[k][j];
                        tt = V[k][i];
                        V[k][i] =  co*tt+si*V[k][j];
                        V[k][j] = -si*tt+co*V[k][j];
                    }
                    // Modify diagonal terms
                    A[i][i] =  co*A[i][i]+si*A[j][i];
                    A[j][j] = -si*A[i][j]+co*A[j][j];
                    A[i][j] = 0.0;
                    // Make 'A' matrix symmetrical
                    for (k=0;k<n;k++) {
                        A[i][k] = A[k][i];
                        A[j][k] = A[k][j];
                    }
                    // A[i][j] made zero by rotation
                }
            }
        }
    } while(Math.abs(ssum)/sum>JXG.Math.eps);
    return [A,V];
};

/**
 * Calculates the integral of function f over interval using Newton-Cotes-algorithm.
 * @param interval The integration interval, e.g. [0, 3]. 
 * @param f A function which takes one argument of type number and returns a number.
 * @return Integral value of f over interval interval
 * @throws {Exception} If {@link JXG.Math.Numerics.number_of_nodes} doesn't match
 * {@link JXG.Math.Numerics.integration_type} an exception is thrown. If you want to use
 * simpson rule respectively milne rule {@link JXG.Math.Numerics.number_of_nodes} must be dividable by
 * 2 respectively 4.
 * @example
 * function f(x) {
 *   return x*x;
 * }
 * 
 * // calculates integral of <tt>f</tt> from 0 to 2.
 * var area1 = JXG.Math.Numerics.NewtonCotes([0, 2], f);
 * 
 * // the same with anonymous function
 * var area2 = JXG.Math.Numerics.NewtonCotes([0, 2], function (x) { return x*x; });
 */
JXG.Math.Numerics.NewtonCotes = function(/** array */ interval, /** function */ f) /** float */ {
    var integral_value = 0.0,
        step_size = (interval[1] - interval[0]) / this.number_of_nodes,
        evaluation_point, i, number_of_intervals;

    switch(this.integration_type) {
        case JXG.INT_TRAPEZ:
            integral_value = (f(interval[0]) + f(interval[1])) * 0.5;
    
            evaluation_point = interval[0];
            for (i = 0; i < this.number_of_nodes - 1; i++)
            {
                evaluation_point += step_size;
                integral_value   += f(evaluation_point);
            }
            integral_value *= step_size;

            break;
        case JXG.INT_SIMPSON:
            if (this.number_of_nodes%2 > 0) {
                throw new Error("JSXGraph:  INT_SIMPSON requires JXG.Math.Numerics.number_of_nodes dividable by 2.");
            }
            number_of_intervals = this.number_of_nodes / 2.0;
            integral_value = f(interval[0]) + f(interval[1]);
            evaluation_point = interval[0];
            for (i = 0; i < number_of_intervals - 1; i++)
            {
                evaluation_point += 2.0 * step_size;
                integral_value   += 2.0 * f(evaluation_point);
            }
            evaluation_point = interval[0] - step_size;
            for (i = 0; i < number_of_intervals; i++)
            {
                evaluation_point += 2.0 * step_size;
                integral_value   += 4.0 * f(evaluation_point);
            }
            integral_value *= step_size / 3.0;
            break;
        default:
            if (this.number_of_nodes%4 > 0) {
                throw new Error("JSXGraph: Error in INT_MILNE: JXG.Math.Numerics.number_of_nodes must be a multiple of 4");
            }
            number_of_intervals = this.number_of_nodes * 0.25;
            integral_value = 7.0 * (f(interval[0]) + f(interval[1]));
            evaluation_point = interval[0];
            for (i = 0; i < number_of_intervals - 1; i++)
            {
                evaluation_point += 4.0 * step_size;
                integral_value   += 14.0 * f(evaluation_point);
            }
            evaluation_point = interval[0] - 3.0 * step_size;
            for (i = 0; i < number_of_intervals; i++)
            {
                evaluation_point += 4.0 * step_size;
                integral_value   += 32.0 * (f(evaluation_point) + f(evaluation_point + 2 * step_size));
            }
            evaluation_point = interval[0] - 2.0 * step_size;
            for (i = 0; i < number_of_intervals; i++)
            {
                evaluation_point += 4.0 * step_size;
                integral_value   += 12.0 * f(evaluation_point);
            }
            integral_value *= 2.0 * step_size / 45.0; /* todo */
    }
    return integral_value;
};

/**
 * Calculates second derivatives at the knots.
 * @param x x values of knots
 * @param y y values of knots
 * @return Second derivatives of interpolated function at the knots.
 */
JXG.Math.Numerics.splineDef = function(/** JXG.Math.Vector */ x, /** JXG.Math.Vector */ y) /** JXG.Math.Vector */ {
    var n = x.length,
        pair, i, diag, z, l,
        data = new Array(),
        dx = [], 
        delta = [],
        F;
        
    if (x.length != y.length)
        throw new Error("JSXGraph: Error in JXG.Math.Numerics.splineDef: Input vector dimensions do not match.");
    
    for (i=0; i<n; i++) {
        pair = {X: x[i], Y: y[i]};
        data.push(pair);
    }
    data.sort(function (a,b) { return a.X - b.X; });
    for (i=0; i<n; i++) {
        x[i] = data[i].X;
        y[i] = data[i].Y;
    }
    
    for (i=0; i<n-1; i++) {
        dx.push(x[i+1] - x[i]);
    }
    for (i=0; i<n-2; i++) {
        delta.push(6 * (y[i+2] - y[i+1])/(dx[i+1]) - 6 * (y[i+1] - y[i])/(dx[i]));
    }

    // ForwardSolve
    diag = new Array();
    z = new Array();
    diag.push(2*(dx[0] + dx[1]));
    z.push(delta[0]);

    for (i=0; i<n-3; i++) {
        l = dx[i+1]/diag[i];
        diag.push(2 * (dx[i+1] + dx[i+2]) - l*dx[i+1]);
        z.push(delta[i+1] - l*z[i]);
    }

    // BackwardSolve
    F = new Array();
    F[n-3] = z[n-3]/diag[n-3];
    for(i=n-4; i>=0; i--) {
        F[i] = (z[i] - (dx[i+1]*F[i+1]))/diag[i];
    }
    
    // Generate f''-Vector
    for(i=n-3; i>=0; i--)
        F[i+1] = F[i];
    
    // natural cubic spline
    F[0] = 0;
    F[n-1] = 0;
    return F;
    //return new JXG.Math.Vector(F);
};

/**
 * Evaluate points on spline.
 * @param {float,Array} x0 A single float value or an array of values to evaluate
 * @param {JXG.Math.Vector} x x values of knots
 * @param {JXG.Math.Vector} y y values of knots
 * @param {JXG.Math.Vector} F Second derivatives at knots, calculated by #splineDef
 * @see splineDef
 * @type float,Array
 * @return A single value
 */
JXG.Math.Numerics.splineEval = function(x0, x, y, F) {
    var n = x.length,
        l = 1,
        asArray = false,
        y0, i, j, a, b, c, d, x_;

    if (n != y.length)
        throw new Error("JSXGraph: Error in JXG.Math.Numerics.splineEval: Defining vector dimensions do not match.");
    
    // number of points to be evaluated
    if(JXG.isArray(x0)) {
        l = x0.length;
        asArray = true;
    } else
        x0 = [x0];
    
    y0 = new Array();
    
    for (i=0; i<l; i++) {
        // is x0 in defining interval?
        if( (x0[i] < x[0]) || (x[i] > x[n-1]))
            return 'NaN';
//            throw new Error("JSXGraph: Error in JXG.Math.Numerics.splineEval: Evaluation point outside spline interval.");
        
        // determine part of spline in which x0 lies
        j;
        for (j=1; j<n; j++) {
            if (x0[i] <= x[j])
                break;
        }
        j--;
        
        // we're now in the j-th partial interval, i.e. x[j] < x0[i] <= x[j+1];
        // determine the coefficients of the polynomial in this interval
        a = y[j];
        b = (y[j+1]-y[j])/(x[j+1]-x[j]) - (x[j+1]-x[j])/6 * (F[j+1]+2*F[j]);
        c = F[j]/2;
        d = (F[j+1]-F[j])/(6*(x[j+1]-x[j]));
        // evaluate x0[i]
        x_ = x0[i]-x[j];
        //y0.push(a + b*x_ + c*x_*x_ + d*x_*x_*x_);
        y0.push(a + (b + (c+ d*x_)*x_)*x_);
    }

    if (asArray)
        return y0;
    else
        return y0[0];
};

/**
  * Generate a string containing the function term of a polynomial.
  * @param {Array} coeffs Coefficients of the polynomial. The position i belongs to x^i.
  * @param {int} deg Degree of the polynomial
  * @param {String} varname Name of the variable (usually 'x')
  * @param {int} prec Precision
  * @return {String} String containg the function term of the polynomial.
  */
JXG.Math.Numerics.generatePolynomialTerm = function(coeffs,deg,varname,prec) {
    var t = '', i;
    for (i=deg;i>=0;i--) {
        t += '('+coeffs[i].toPrecision(prec)+')';
        if (i>1) { t+='*'+varname+'<sup>'+i+'</sup> + '; }
        else if (i==1) { t+='*'+varname+' + '; }
    }
    return t;
}

/**
 * Computes the polynomial through a given set of coordinates in Lagrange form.
 * Returns the Lagrange polynomials, see
 * Jean-Paul Berrut, Lloyd N. Trefethen: Barycentric Lagrange Interpolation,
 * SIAM Review, Vol 46, No 3, (2004) 501-517.
 * @param {array} p Array of JXG.Points
 * @type {function}
 * @return {function} A function of one parameter which returns the value of the polynomial,
 * whose graph runs through the given points.
 */
JXG.Math.Numerics.lagrangePolynomial = function(p) {  
    var w = [];
    var term = '';
    var fct = function(x,suspendedUpdate)  {
        var i, k, len, xi, s,
            num = 0, denom = 0;
        
        len = p.length;
        if (!suspendedUpdate) {
            for (i=0;i<len;i++) {
                w[i] = 1.0;
                xi = p[i].X();
                for (k=0;k<len;k++) if (k!=i) {
                    w[i] *= (xi-p[k].X());
                }
                w[i] = 1/w[i];
             }

            M = [];
            for (j=0;j<len;j++) {
                M.push([1]);
            }
                /* // Function term not yet
                for (i=1;i<len;i++) {
                    for (j=0;j<len;j++) {
                        M[j][i] = M[j][i-1]*datax[j];      // input data
                    }
                }
                y = curve.dataY;                           // input data
                MT = JXG.Math.Matrix.transpose(M);

                B = JXG.Math.matMatMult(MT,M);
                c = JXG.Math.matVecMult(MT,y);
                container = JXG.Math.Numerics.Gauss(B, c);
                
                for (i=degree, term='';i>=0;i--) {
                    term += '('+container[i].toPrecision(3)+')';
                    if (i>1) { term+='*x<sup>'+i+'</sup> + '; }
                    else if (i==1) { term+='*x + '; }
                }
                */
        }
        
        for (i=0;i<len;i++) {
            xi = p[i].X();
            if (x==xi) { 
                return p[i].Y(); 
            } else {
                s = w[i]/(x-xi)
                denom += s;
                num += s*p[i].Y();
            }
        }
        return num/denom;
    }
    fct.getTerm = function() {
        return term;
    };
    
    return fct;

/*
    return function(x) {
        var i,k,t,
            len = p.length,
            y = 0.0,
            xc = [];
        
        for (i=0;i<len;i++) {
            xc[i] = p[i].X();
        }
        for (i=0;i<len;i++) {
            t = p[i].Y();
            for (k=0;k<len;k++) if (k!=i) {
                t *= (x-xc[k])/(xc[i]-xc[k]);
            }
            y += t;
        }
        return y;
    };
*/    
};

/**
 * Returns the Lagrange polynomials for curves with equidistant nodes, see
 * Jean-Paul Berrut, Lloyd N. Trefethen: Barycentric Lagrange Interpolation,
 * SIAM Review, Vol 46, No 3, (2004) 501-517.
 * The graph of the parametric curve [f(t),g(t)] runs through the given points.
 * @param {Array} p Artray of JXG.Points
 * @type {Array function, function value, value]}
 * @return {array} [f(t),g(t),0,p.length-1],
 */
JXG.Math.Numerics.neville = function(p) {
    var w = [];
    
    var xfct = function(t, suspendedUpdate) {
        var i, d, L, s, 
            bin = JXG.Math.binomial,
            len = p.length,
            len1 = len - 1,
            num = 0.0, 
            denom = 0.0;
            
        if (!suspendedUpdate) {
            s = 1;
            for (i=0;i<len;i++) {
                w[i] = bin(len1,i)*s;
                s *= (-1);
            }
        }

        d = t;
        for (i=0;i<len;i++) {
            if (d==0) {
                return p[i].X();
            } else {
                s = w[i]/d;
                d--;
                num   += p[i].X()*s;
                denom += s;
            }
        }
        return num/denom;
    }
    var yfct = function(t, suspendedUpdate) {
        var i, d, L, s, 
            bin = JXG.Math.binomial,
            len = p.length,
            len1 = len - 1,
            num = 0.0, 
            denom = 0.0;
            
        if (!suspendedUpdate) {
            //L = JXG.Math.binomial(len-1,i)*((i%2==0)?1:(-1))/d;
            s = 1;
            for (i=0;i<len;i++) {
                w[i] = bin(len1,i)*s;
                s *= (-1);
            }
        }

        d = t;
        for (i=0;i<len;i++) {
            if (d==0) {
                return p[i].Y();
            } else {
                s = w[i]/d;
                d--;
                num   += p[i].Y()*s;
                denom += s;
            }
        }
        return num/denom;
    }
    return [xfct, yfct, 0, function(){ return p.length-1;}];
};

/**
 * Computes the regression polynomial of a given degree through a given set of coordinates.
 * Returns the regression polynomial function.
 * @param degree number, function or slider.
 * Either
 * @param dataX array containing the x-coordinates of the data set
 * @param dataY array containing the y-coordinates of the data set, 
 * or
 * @param data array consisting of JXG.Points.
 * @type {function}
 * @return {function} A function of one parameter which returns the value of the regression polynomial of the given degree.
 * It possesses the method getTerm() which returns the string containing the function term of the polynomial.
 */
JXG.Math.Numerics.regressionPolynomial = function(degree, dataX, dataY) { 
    var coeffs = [],
        dbg_count = 0,
        deg, dX, dY,
        inputType,
        term = '';
    
    if (JXG.isPoint(degree) && typeof degree.Value == 'function') {  // Slider
        deg = function(){return degree.Value();};
    } else if (JXG.isFunction(degree)) {
        deg = degree;
    } else if (JXG.isNumber(degree)) {
        deg = function(){return degree;};
    } else {
        throw new Error("JSXGraph: Can't create regressionPolynomial from degree of type'" + (typeof degree) + "'.");
    }
    
    if (arguments.length==3 && JXG.isArray(dataX) && JXG.isArray(dataY)) {              // Parameters degree, dataX, dataY
//        dX = dataX;
//        dY = dataY;
// is done later
        inputType = 0;
    } else if ( arguments.length==2 && JXG.isArray(dataX) && JXG.isPoint(dataX[0]) ) {  // Parameters degree, point array
        inputType = 1;
    } else {
        throw new Error("JSXGraph: Can't create regressionPolynomial. Wrong parameters.");
    }
    
    var fct = function(x,suspendedUpdate){
            var i, j, M, MT, y, B, c, s,
                d,
                len = dataX.length;                        // input data
                
            d = Math.floor(deg());                      // input data
            if (!suspendedUpdate) {
                if (inputType==1) {  // point list as input 
                    dX = [];
                    dY = [];
                    for (i=0;i<len;i++) {
                        dX[i] = dataX[i].X();
                        dY[i] = dataX[i].Y();
                    }
                }
                
                if (inputType==0) {  // check for functions
                    dX = [];
                    dY = [];
                    for(i=0;i<len;i++) {
                        if(JXG.isFunction(dataX[i]))
                            dX.push(dataX[i]());
                        else
                            dX.push(dataX[i]);
                        if(JXG.isFunction(dataY[i]))
                            dY.push(dataY[i]());
                        else
                            dY.push(dataY[i]);
                    }
                }

                M = [];
                for (j=0;j<len;j++) {
                    M.push([1]);
                }
                for (i=1;i<=d;i++) {
                    for (j=0;j<len;j++) {
                        M[j][i] = M[j][i-1]*dX[j];      // input data
                    }
                }
                
                y = dY;                                 // input data
                MT = JXG.Math.Matrix.transpose(M);
                B = JXG.Math.matMatMult(MT,M);
                c = JXG.Math.matVecMult(MT,y);
                coeffs = JXG.Math.Numerics.Gauss(B, c);

                term = JXG.Math.Numerics.generatePolynomialTerm(coeffs,d,'x',3);         
            }
            
            // Horner's scheme to evaluate polynomial
            s = coeffs[d];
            for (i=d-1;i>=0;i--) {
                s = (s*x+coeffs[i]);
            }
            return s;
    };
    fct.getTerm = function() {
        return term;
    };
    return fct;
}
    
/**
 * Computes the cubic Bezier curve through a given set of points..
 * @param data array consisting of 3*k+1 JXG.Points.
 * The points at position k with k mod 3 = 0 are the data points,
 * points at position k with k mod 3 = 1 or 2 are the control points.
 * @type {function}
 * @return {function} A function of one parameter which returns the value of the cubic Bezier curve.
 */
JXG.Math.Numerics.bezier = function(points) {
    var len = 0; 

    return [function(t,suspendedUpdate) {
                var z = Math.floor(t)*3,
                    t0 = t % 1,
                    t1 = 1-t0;
                        
                if (!suspendedUpdate) {
                    len = Math.floor(points.length/3);
                }
                        
                if (t<0) { return points[0].X(); }
                if (t>=len) { return points[points.length-1].X(); }
                if (isNaN(t)) { return NaN; }
                return t1*t1*(t1*points[z].X()+3*t0*points[z+1].X())+(3*t1*points[z+2].X()+t0*points[z+3].X())*t0*t0;
            },
            function(t,suspendedUpdate) {
                var z = Math.floor(t)*3,
                    t0 = t % 1,
                    t1 = 1-t0;
                        
                if (!suspendedUpdate) {
                    len = Math.floor(points.length/3);
                }
                        
                if (t<0) { return points[0].Y(); }
                if (t>=len) { return points[points.length-1].Y(); }
                if (isNaN(t)) { return NaN; }
                return t1*t1*(t1*points[z].Y()+3*t0*points[z+1].Y())+(3*t1*points[z+2].Y()+t0*points[z+3].Y())*t0*t0;
            }, 
            0, function() {return Math.floor(points.length/3);}];
};

/**
 * Numerical (symmetric) approximation of derivative.
 * @param {function} f Function in one variable to be differentiated.
 * @param {object} obj Optional object that is treated as "this" in the function body. This is useful, if the function is a 
 *                 method of an object and contains a reference to its parent object via "this".
 * suspendUpdate is piped through, {@link JXG.Curve#updateCurve} and {@link JXG.Curve#hasPoint}.
 * @type {function}
 * @return {function} Derivative function of a given function f.
 */
JXG.Math.Numerics.D = function(/** function */ f, /** object */ obj) /* function */ {
    var h = 0.00001,
        h2 = 1.0/(h*2.0);
    
    if (arguments.length==1 || (arguments.length>1 && typeof arguments[1]=='undefined') ){ 
        return function(x,suspendUpdate){ return (f(x+h,suspendUpdate)-f(x-h,suspendUpdate))*h2; };
    } else {                   // set "this" to "obj" in f 
        return function(x,suspendUpdate){ return (f.apply(obj,[x+h,suspendUpdate])-f.apply(obj,[x-h,suspendUpdate]))*h2; };
    }
};

/**
 * Integral of function f over interval.
 * @deprecated Use {@link JXG.Math.Numerics.NewtonCotes} instead.
 * @see JXG.Math.Numerics.NewtonCotes
 */
JXG.Math.Numerics.I = function(/** array */ interval, /** function */ f) {
    return JXG.Math.Numerics.NewtonCotes(interval, f);
};

/**
 * Newton's method to find roots of a funtion in one variable.
 * @param {function} f We search for a solution of f(x)=0.
 * @param {float} x initial guess for the root, i.e. staring value.
 * @param {object} obj optional object that is treated as "this" in the function body. This is useful, if the function is a 
 *                 method of an object and contains a reference to its parent object via "this".
 * @return {float} root of the function f.
 */
JXG.Math.Numerics.newton = function(/** function */ f, /** number */ x, /** object */ obj) /** number */ {
    var i = 0,
        h = 0.000001,
        newf = f.apply(obj,[x]), // set "this" to "obj" in f 
        df;
        
    while (i<50 && Math.abs(newf)>h) {
        df = this.D(f,obj)(x);
        if (Math.abs(df)>h) {
            x -= newf/df;
        } else {
            x += (Math.random()*0.2-1.0);
        }
        newf = f.apply(obj,[x]);
        i++;
    }
    return x;
};

/**
 * Abstract method to find roots of univariate functions.
 * @param {function} f We search for a solution of f(x)=0.
 * @param {float} x initial guess for the root, i.e. staring value.
 * @param {object} obj optional object that is treated as "this" in the function body. This is useful, if the function is a 
 *                 method of an object and contains a reference to its parent object via "this".
 * @return {float} root of the function f.
 */
JXG.Math.Numerics.root = function(/** function */ f, /** number */ x, /** object */ obj) /** number */ {
    return this.newton(f,x,obj);
};

/**
 * Hlper function to create curve which displays Riemann sums.
 * Compute coordinates for the rectangles showing the Riemann sum.
 * @param {function} f Function f, whose integral is approximated by the Riemann sum.
 * @param {int} n number of rectangles.
 * @param {String} type Type of approximation. Possible values are: 'left', 'right', 'middle', 'lower', 'upper', or 'trapezodial'.
 * @param {float} start Left border of the approximation interval
 * @param {float} end Right border of the approximation interval
 * @return {array} An array of two arrays containing the x and y coordinates for the rectangles showing the Riemann sum. This array may be used as
 *                 parent array of a JXG.Curve.
 */
JXG.Math.Numerics.riemann = function(/** function */ f, /** type */ n,  /** type */ type,  /** type */ start,  /** type */ end)  /** array */ {
    var xarr,yarr,i,delta,j,x,y,x1,delta1,y1;
    
    xarr = [];
    yarr = [];
    j = 0;
    x = start;
    n = Math.floor(n);
    xarr[j] = x; yarr[j] = 0.0;
    
    if (n>0) {
        delta = (end-start)/n;
        delta1 = delta*0.01; // for 'lower' and 'upper'
        
        for (i=0;i<n;i++) {
            if (type=='right') {
                y = f(x+delta);
            } else if (type=='middle') {
                y = f(x+delta*0.5);
            } else if ((type=='left') || (type=='trapezodial')) {
                y = f(x);
            } else if (type=='lower') {
                y = f(x);
                for (x1=x+delta1;x1<=x+delta;x1+=delta1) {
                    y1 = f(x1);
                    if (y1<y) { y = y1; };
                }
            } else { // (type=='upper')
                y = f(x);
                for (x1=x+delta1;x1<=x+delta;x1+=delta1) {
                    y1 = f(x1);
                    if (y1>y) { y = y1; };
                }
            }
            
            j++;
            xarr[j] = x; yarr[j] = y;
            j++; x+=delta;
            if (type=='trapezodial') {
                y = f(x);
            }
            xarr[j] = x; yarr[j] = y;
            j++;
            xarr[j] = x; yarr[j] = 0.0;
         }
    }
    return [xarr,yarr];
};

/**
 * Approximate the integral by Riemann sums.
 * Compute the area described by the riemann sum rectangles.
 * @param {function} f Function f, whose integral is approximated by the Riemann sum.
 * @param {int} n number of rectangles.
 * @param {String} type Type of approximation. Possible values are: 'left', 'right', 'middle', 'lower', 'upper', or 'trapezodial'.
 * @param {float} start Left border of the approximation interval
 * @param {float} end Right border of the approximation interval
 * @return {float} The sum of the areas of the rectangles.
 */
JXG.Math.Numerics.riemannsum = function(/** function */ f, /** type */ n,  /** type */ type,  /** type */ start,  /** type */ end)  /** number */ {
    var sum,i,delta,x,y,x1,delta1,y1;
    
    sum = 0.0;
    x = start;
    n = Math.floor(n);
    if (n>0) {
        delta = (end-start)/n;
        delta1 = delta*0.01; // for 'lower' and 'upper'
        for (i=0;i<n;i++) {
            if (type=='right') {
                y = f(x+delta);
            } else if (type=='middle') {
                y = f(x+delta*0.5);
            } else if (type=='trapezodial') {
                y = 0.5*(f(x+delta)+f(x));
            } else if (type=='left') { 
                y = f(x);
            } else if (type=='lower') {
                y = f(x);
                for (x1=x+delta1;x1<=x+delta;x1+=delta1) {
                    y1 = f(x1);
                    if (y1<y) { y = y1; };
                }
            } else { // (type=='upper')
                y = f(x);
                for (x1=x+delta1;x1<=x+delta;x1+=delta1) {
                    y1 = f(x1);
                    if (y1>y) { y = y1; };
                }
            }
            sum += delta*y;
            x += delta;
         }
    }
    return sum;
};

/**
 * Object for storing butcher tableaus for Runge-Kutta-methods.
 * @class
 * @description
 * @see JXG.Math.Numerics.rungeKutta
 */
JXG.Math.Numerics.Butcher = function () {
    /**
     * Order of Runge-Kutta-method.
     * @type number
     */
    this.s = 0;

    /**
     * 2-dimensional array containing the butcher tableau matrix.
     * See <a href="http://en.wikipedia.org/wiki/Runge-Kutta_methods">http://en.wikipedia.org/wiki/Runge-Kutta_methods</a>.
     * @type array
     */
    this.A = [];

    /**
     * Array containing the coefficients below the butcher tableau matrix.
     * See <a href="http://en.wikipedia.org/wiki/Runge-Kutta_methods">http://en.wikipedia.org/wiki/Runge-Kutta_methods</a>.
     * @type array
     */
    this.b = [];

    /**
     * Array containing the coefficients to the left of the butcher tableau matrix.
     * See <a href="http://en.wikipedia.org/wiki/Runge-Kutta_methods">http://en.wikipedia.org/wiki/Runge-Kutta_methods</a>.
     * @type array
     */
    this.c = [];
};

/**
 * Predefined butcher tableaus for the common Runge-Kutta method (fourth order), Heun method (second order), and Euler method (first order).
 * @namespace
 */
JXG.Math.Numerics.predefinedButcher = {};

/**
 * Butcher tableau for common fourth order Runge-Kutta method.
 * @type JXG.Math.Numerics.Butcher
 */
JXG.Math.Numerics.predefinedButcher.RK4 = {
    s: 4,
    A: [[ 0,  0,  0, 0],
        [0.5, 0,  0, 0],
        [ 0, 0.5, 0, 0],
        [ 0,  0,  1, 0]],
    b: [1./6., 1./3., 1./3., 1./6.],
    c: [0, 0.5, 0.5, 1]
};

/**
 * Butcher tableau for heun method.
 * @type JXG.Math.Numerics.Butcher
 */
JXG.Math.Numerics.predefinedButcher.Heun = {
    s: 2,
    A: [[0, 0], [1, 0]],
    b: [0.5, 0.5],
    c: [0, 1]
};

/**
 * Butcher tableau for euler method.
 * @type JXG.Math.Numerics.Butcher
 */
JXG.Math.Numerics.predefinedButcher.Euler = {
    s: 1,
    A: [[0]],
    b: [1],
    c: [0]
};

/**
 * Solve initial value problems numerically using Runge-Kutta-methods.
 * See {@link http://en.wikipedia.org/wiki/Runge-Kutta_methods} for more information on the algorithm.
 * @param butcher Butcher tableau describing the Runge-Kutta method to use.
 * @param x0 Initial value vector. If the problem is of one-dimensional, the initial value also has to be given in an array.
 * @param I Interval on which to integrate.
 * @param N Number of evaluation points.
 * @param f Function describing the right hand side of the first order ordinary differential equation, i.e. if the ode
 * is given by the equation <pre>dx/dt = f(t, x(t)).</pre> So f has to take two parameters, a number <tt>t</tt> and a
 * vector <tt>x</tt>, and has to return a vector of the same dimension as <tt>x</tt> has.
 * @return An array of vectors describing the solution of the ode on the given interval I.
 * @example
 * // A very simple autonomous system dx(t)/dt = x(t);
 * function f(t, x) {
 *     return x;
 * }
 * 
 * // We want to use the method of heun.
 * var method = JXG.Math.Numerics.predefinedButcher.Heun;
 * // Solve it with initial value x(0) = 1 on the interval [0, 2]
 * // with 20 evaluation points.
 * var data = JXG.Math.Numerics.rungeKutta(method, [1], [0, 2], 20, f);
 * 
 * // Prepare data for plotting the solution of the ode using a curve. 
 * var dataX = [];
 * var dataY = [];
 * var h = 0.1;        // (I[1] - I[0])/N  = (2-0)/20
 * for(var i=0; i&lt;data.length; i++) {
 *     dataX[i] = i*h;
 *     dataY[i] = data[i][0];
 * }
 * var g = board.create('curve', [dataX, dataY], {strokeWidth:'2px'});
 * </pre><div id="d2432d04-4ef7-4159-a90b-a2eb8d38c4f6" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 * var board = JXG.JSXGraph.initBoard('d2432d04-4ef7-4159-a90b-a2eb8d38c4f6', {boundingbox: [-1, 5, 5, -1], axis: true, showcopyright: false, shownavigation: false});
 * function f(t, x) {
 *     // we have to copy the value.
 *     // return x; would just return the reference.
 *     return [x[0]];
 * }
 * var data = JXG.Math.Numerics.rungeKutta(JXG.Math.Numerics.predefinedButcher.Heun, [1], [0, 2], 20, f);
 * var dataX = [];
 * var dataY = [];
 * var h = 0.1;
 * for(var i=0; i<data.length; i++) {
 *     dataX[i] = i*h;
 *     dataY[i] = data[i][0];
 * }
 * var g = board.create('curve', [dataX, dataY], {strokeColor:'red', strokeWidth:'2px'});
 * </script><pre>
 */
JXG.Math.Numerics.rungeKutta = function(/** JXG.Math.Numerics.Butcher */ butcher, /** array */ x0,
										/** array */ I, /** number */ N, /** function */ f) /** array */ {

	// TODO error/parameter check:
    // N not too big (warn or give up?) OR adaptive stepsize.

    var x = [],
        y = [],
        h = (I[1]-I[0])/N,
        t = I[0],
        e, i, j, 
        k, l,
        dim = x0.length,
        s = butcher.s,
        numberOfResultPoints = 1000,
        quotient = N/numberOfResultPoints,
        result = [],
        r = 0;
        
    // don't change x0, so copy it
    for (e=0; e<dim; e++)
        x[e] = x0[e];
    for (i=0; i<N; i++) {
        // Optimization doesn't work for ODEs plotted using time
//        if((i % quotient == 0) || (i == N-1)) {
            result[r] = [];
            for (e=0; e<dim; e++)
                result[r][e] = x[e];
            r++;
//        }
        // init k
        k = [];
        for(j=0; j<s; j++) {
            // init y = 0
            for (e=0; e<dim; e++)
                y[e] = 0.;

            // Calculate linear combination of former k's and save it in y
            for (l=0; l<j; l++) {
                for (e=0; e<dim; e++) {
                    y[e] += (butcher.A[j][l])*h*k[l][e];
                }
            }

            // add x(t) to y
            for(e=0; e<dim; e++) {
                y[e] += x[e];
            }

            // calculate new k and add it to the k matrix
            k.push(f(t+butcher.c[j]*h, y));
        }

        // init y = 0
        for (e=0; e<dim; e++)
            y[e] = 0.;

        for (l=0; l<s; l++) {
            for (e=0; e<dim; e++)
                y[e] += butcher.b[l]*k[l][e];
        }

        for (e=0; e<dim; e++) {
            x[e] = x[e] + h*y[e];
        }

        t += h;
    }

    return result;
};

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
  * Functions for mathematical statistics
  * Most functions are R-like:
  * For example prod(a1,a2) computes an array c such that
  * for (i=0;i<a1.length;i++) c[i] = a1[i]*a2[i];
  *
 **/
JXG.Math.Statistics = {};

JXG.Math.Statistics.sum = function(arr) {
    var i, len, res = 0;
    
    for(i=0, len=arr.length; i<len; i++) { 
        res += arr[i];
    } 
    return res;
};

JXG.Math.Statistics.prod = function(arr) {
    var i, len, res = 1;
    
    for(i=0, len=arr.length; i<len; i++) { 
        res *= arr[i];
    } 
    return res;
};

JXG.Math.Statistics.mean = function(arr) {
    if (arr.length>0) {
        return this.sum(arr)/arr.length;
    } else {
        return 0.0;
    }
};

JXG.Math.Statistics.median = function(arr) {
    var tmp, len;
    
    if (arr.length>0) {
        tmp = arr.clone();
        tmp.sort(function(a,b){return a-b;});
        len = tmp.length;
        if (len%2==1) {
            return tmp[parseInt(len*0.5)];
        } else{
            return (tmp[len*0.5-1]+tmp[len*0.5])*0.5;
        }
    } else {
        return 0.0;
    }
};

/**
 * bias-corrected sample variance
 */
JXG.Math.Statistics.variance = function(arr) {
    var m, res, i, len;
    
    if (arr.length>1) {
        m = this.mean(arr);
        res = 0;
        for(i=0, len=arr.length; i<len; i++) { 
            res += (arr[i]-m)*(arr[i]-m);
        } 
        return res/(arr.length-1);
    } else {
        return 0.0;
    }
};

JXG.Math.Statistics.sd = function(arr) {
    return Math.sqrt(this.variance(arr));
};

JXG.Math.Statistics.weightedMean = function(arr,w) {
    if (arr.length!=w.length) { return; }
    if (arr.length>0) {
        return this.mean(this.multiply(arr,w));
    } else {
        return 0.0;
    }
};

JXG.Math.Statistics.max = function(arr) {
    var res, i, len;
    
    if (arr.length==0) { return NaN; }
    res = arr[0];
    for(i=1, len=arr.length; i<len; i++) { 
        res = (arr[i]>res)?(arr[i]):res;
    } 
    return res;
};

JXG.Math.Statistics.min = function(arr) {
    var res, i, len;

    if (arr.length==0) { return NaN; }
    res = arr[0];
    for(i=1, len=arr.length; i<len; i++) { 
        res = (arr[i]<res)?(arr[i]):res;
    } 
    return res;
};

/**
 * R-style functions
 */
JXG.Math.Statistics.range = function(arr) {
    return [this.min(arr),this.max(arr)];
};

JXG.Math.Statistics.diff = function(arr) { // ?????
    return arr;
};

JXG.Math.Statistics.min = function(arr) {
    var res, i, len;
    
    if (arr.length==0) { return NaN; }
    res = arr[0];
    for(i=1, len=arr.length; i<len; i++) { 
        res = (arr[i]<res)?(arr[i]):res;
    } 
    return res;
};

JXG.Math.Statistics.abs = function(arr) {  // This can be generalized with Prototype.js and should be done for all Math. methods
    var i, len, res = [];
    if (typeof JXG.isArray(arr1)) {
        for (i=0, len=arr.length;i<len;i++) { res[i] = Math.abs(arr[i]); }
    } else if (typeof arr=='number') {
        return Math.abs(arr);
    } else {
        res = null;
    }
    return res;
};

JXG.Math.Statistics.add = function(arr1,arr2) {
    var i, len, res = [];
    
    if (typeof JXG.isArray(arr1) && typeof arr2=='number') {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1[i]+arr2; }
    } else if (typeof arr1=='number' && typeof JXG.isArray(arr2)) {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1+arr2[i]; }
    } else if (typeof JXG.isArray(arr1) && typeof JXG.isArray(arr2)) {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1[i]+arr2[i]; }
    } else if (typeof arr1=='number' && typeof arr2=='number') {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1+arr2; }
    } else {
        res = null;
    }
    return res;
};

JXG.Math.Statistics.divide = function(arr1,arr2) {
    var i, len, res = [];
    
    if (typeof JXG.isArray(arr1) && typeof arr2=='number') {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1[i]/arr2; }
    } else if (typeof arr1=='number' && typeof JXG.isArray(arr2)) {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1/arr2[i]; }
    } else if (typeof JXG.isArray(arr1) && typeof JXG.isArray(arr2)) {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1[i]/arr2[i]; }
    } else if (typeof arr1=='number' && typeof arr2=='number') {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1/arr2; }
    } else {
        res = null;
    }
    return res;
};

JXG.Math.Statistics.mod = function(arr1,arr2) {
    var i, len, res = [];
    
    if (typeof JXG.isArray(arr1) && typeof arr2=='number') {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1[i]%arr2; }
    } else if (typeof arr1=='number' && typeof JXG.isArray(arr2)) {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1%arr2[i]; }
    } else if (typeof JXG.isArray(arr1) && typeof JXG.isArray(arr2)) {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1[i]%arr2[i]; }
    } else if (typeof arr1=='number' && typeof arr2=='number') {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1%arr2; }
    } else {
        res = null;
    }
    return res;
};

JXG.Math.Statistics.multiply = function(arr1,arr2) {
    var i, len, res = [];
    
    if (typeof JXG.isArray(arr1) && typeof arr2=='number') {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1[i]*arr2; }
    } else if (typeof arr1=='number' && typeof JXG.isArray(arr2)) {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1*arr2[i]; }
    } else if (typeof JXG.isArray(arr1) && typeof JXG.isArray(arr2)) {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1[i]*arr2[i]; }
    } else if (typeof arr1=='number' && typeof arr2=='number') {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1*arr2; }
    } else {
        res = null;
    }
    return res;
};

JXG.Math.Statistics.subtract = function(arr1,arr2) {
    var i, len, res = [];
    
    if (typeof JXG.isArray(arr1) && typeof arr2=='number') {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1[i]-arr2; }
    } else if (typeof arr1=='number' && typeof JXG.isArray(arr2)) {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1-arr2[i]; }
    } else if (typeof JXG.isArray(arr1) && typeof JXG.isArray(arr2)) {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1[i]-arr2[i]; }
    } else if (typeof arr1=='number' && typeof arr2=='number') {
        for (i=0, len=Math.min(arr1.length,arr2.length);i<len;i++) { res[i] = arr1-arr2; }
    } else {
        res = null;
    }
    return res;
};

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph. If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview In this file the namespace Math.Symbolic is defined, which holds methods
 * and algorithms for symbolic computations.
 * @author graphjs
 */
/**
 * Math.Symbolic
 */
JXG.Math.Symbolic = {};

/**
 * Generates symbolic coordinates for the part of a construction including all the elements from that
 * a specific element depends of. These coordinates will be stored in GeometryElement.symbolic.
 * @param {JXG.Board} board The board that's element get some symbolic coordinates.
 * @param {JXG.GeometryElement} element All ancestor of this element get symbolic coordinates.
 * @param {String} variable Name for the coordinates, e.g. x or u.
 * @param {String} append Method for how to append the number of the coordinates. Possible values are
 *                        'underscore' (e.g. x_2), 'none' (e.g. x2), 'brace' (e.g. x[2]).
 * @type int
 * @return Number of coordinates given.
 */
JXG.Math.Symbolic.generateSymbolicCoordinatesPartial = function(board, element, variable, append) {
    function makeCoords(num) {
        if (append == 'underscore')
            return '' + variable + '_{' + num + '}';
        else if (append == 'brace')
            return '' + variable + '[' + num + ']';
        else
            return '' + variable + '' + num;
    }

    var list = element.ancestors;
    var count = 0;
    var t_num;
    for(var t in list) {
        t_num = 0;
        if (JXG.isPoint(list[t])) {
            for(var k in list[t].ancestors) {
                t_num++;
            }
            if(t_num == 0) {
                list[t].symbolic.x = list[t].coords.usrCoords[1];
                list[t].symbolic.y = list[t].coords.usrCoords[2];
            } else {
                count++;
                list[t].symbolic.x = makeCoords(count);
                count++;
                list[t].symbolic.y = makeCoords(count);
            }

        }
    }

    if(JXG.isPoint(element)) {
        element.symbolic.x = 'x';
        element.symbolic.y = 'y';
    }

    return count;
};

/**
 * Clears all .symbolic.x and .symbolic.y members on every point of a given board.
 * @param {JXG.Board} board The board that's points get cleared their symbolic coordinates.
 */
JXG.Math.Symbolic.clearSymbolicCoordinates = function(board) {
    for(var t in board.objects) {
        if (JXG.isPoint(board.objects[t])) {
            board.objects[t].symbolic.x = '';
            board.objects[t].symbolic.y = '';
        }
    }
};

/**
 * Generates polynomials for the part of a construction including all the points from that
 * a specific element depends of.
 * @param {JXG.Board} board The board that's points polynomials will be generated.
 * @param {JXG.GeometryElement} element All points in the set of ancestors of this element are used to generate the set of polynomials.
 * @type Array
 * @return Array of polynomials as strings.
 */
JXG.Math.Symbolic.generatePolynomials = function(board, element, generateCoords) {
    if(generateCoords)
        this.generateSymbolicCoordinatesPartial(board, element, 'u', 'brace');

    var list = element.ancestors,
        number_of_ancestors,
        pgs = [],
        result = [],
        t, k, i;
    list[element.id] = element;

    for(t in list) {
        number_of_ancestors = 0;
        pgs = [];
        if (JXG.isPoint(list[t])) {
            for(k in list[t].ancestors) {
                number_of_ancestors++;
            }
            if(number_of_ancestors > 0) {
                pgs = list[t].generatePolynomial();
                for(i=0; i<pgs.length; i++)
                    result.push(pgs[i]);
            }
        }
    }

    if(generateCoords)
        this.clearSymbolicCoordinates(board);

    return result;
};

/**
 * Calculate geometric locus of a point given on a board. Invokes python script on server.
 * @param {JXG.Board} board The board on that the point lies.
 * @param {JXG.Point} point The point that will be traced.
 * @param {function} callback A callback function that is called after the server request is finished.
 *    Must take an array of strings as the only parameter.
 * @type Array
 * @return Array of points.
 */
JXG.Math.Symbolic.geometricLocusByGroebnerBase = function(board, point, callback) {
    var numDependent = this.generateSymbolicCoordinatesPartial(board, point, 'u', 'brace'),
        poly = this.generatePolynomials(board, point);
    var polyStr = poly.join(','),
        xsye = new JXG.Coords(JXG.COORDS_BY_USR, [0,0], board),
        xeys = new JXG.Coords(JXG.COORDS_BY_USR, [board.canvasWidth, board.canvasHeight], board),
        fileurl;

    if(typeof JXG.Server.modules.geoloci == 'undefined')
        JXG.Server.loadModule('geoloci')

    if(typeof JXG.Server.modules.geoloci == 'undefined')
        throw new Error("JSXGraph: Unable to load JXG.Server module 'geoloci.py'.");

    this.cbp = function(data) {
        //alert(data.exectime);
        callback(data.datax, data.datay, data.polynomial);
    };

    this.cb = JXG.bind(this.cbp, this);

    JXG.Server.modules.geoloci.lociCoCoA(xsye.usrCoords[1], xeys.usrCoords[1], xeys.usrCoords[2], xsye.usrCoords[2], numDependent, polyStr, this.cb);

    this.clearSymbolicCoordinates(board);
};

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/


/**
 * @fileoverview A class for complex arithmetics JXG.Complex is defined in this
 * file. Also a namespace JXG.C is included to provide instance-independent
 * arithmetic functions.
 * @author graphjs
 */

/**
 * Creates a new complex number.
 * @class This class is for calculating with complex numbers.
 * @param [x=0] Real part of the resulting complex number.
 * @param [y=0] Imaginary part of the resulting complex number.
 * @returns An object representing the complex number <tt>x + iy</tt>.
 */
JXG.Complex = function(/** number */ x, /** number */ y) {
	/**
	 * This property is only to signalize that this object is of type JXG.Complex. Only
	 * used internally to distinguish between normal JavaScript numbers and JXG.Complex numbers.
	 * @type boolean
	 * @default true
	 * @private
	 */
    this.isComplex = true;
    
    if (typeof x == 'undefined') {
        x = 0;
    }
    if (typeof y == 'undefined') {
        y = 0;
    }

    /* is the first argument a complex number? if it is,
     * extract real and imaginary part. */
    if (x.isComplex) {
        y = x.imaginary;
        x = x.real;
    }

    /**
     * Real part of the complex number.
     * @type number
     * @default 0
     */
    this.real = x;

    /**
     * Imaginary part of the complex number.
     * @type number
     * @default 0
     */
    this.imaginary = y;

    /**
     * Absolute value in the polar form of the complex number. Currently unused.
     * @type number
     */
    this.absval = 0;

    /**
     * Angle value in the polar form of the complex number. Currently unused.
     * @type number
     */
    this.angle = 0;
};

/**
 * Converts a complex number into a string.
 * @return Formatted string containing the complex number in human readable form (algebraic form).
 */
JXG.Complex.prototype.toString = function() /** string */{
    return '' + this.real + ' + ' + this.imaginary + 'i';
};

/**
 * Add another complex number to this complex number.
 * @param c A JavaScript number or a JXG.Complex object to be added to the current object.
 */
JXG.Complex.prototype.add = function(/** JXG.Complex,number */ c) /** undefined */ {
    if(typeof c == 'number') {
        this.real += c;
    } else {
        this.real += c.real;
        this.imaginary += c.imaginary;
    }
};

/**
 * Subtract another complex number from this complex number.
 * @param c A JavaScript number or a JXG.Complex object to subtract from the current object.
 */
JXG.Complex.prototype.sub = function(/** JXG.Complex,number */ c) /** undefined */{
    if(typeof c == 'number') {
        this.real -= c;
    } else {
        this.real -= c.real;
        this.imaginary -= c.imaginary;
    }
};

/**
 * Multiply another complex number to this complex number.
 * @param c A JavaScript number or a JXG.Complex object to
 * multiply with the current object.
 */
JXG.Complex.prototype.mult = function(/** JXG.Complex,number */ c) /** undefined */{
    if(typeof c == 'number') {
        this.real *= c;
        this.imaginary *= c;
    } else {
        //  (a+ib)(x+iy) = ax-by + i(xb+ay)
        this.real = this.real*c.real - this.imaginary*c.imaginary;
        this.imaginary = this.real*c.imaginary + this.imaginary*c.real;
    }
};


/**
 * Divide this complex number by the given complex number.
 * @param c A JavaScript number or a JXG.Complex object to
 * divide the current object by.
 */
JXG.Complex.prototype.div = function(/** JXG.Complex,number */ c) /** undefined */{
    var denom;

    if(typeof c == 'number') {
        if(Math.abs(c) < Math.eps) {
            this.real = Infinity;
            this.imaginary = Infinity;
            
            return;
        }
        this.real /= c;
        this.imaginary /= c;
    } else {
        //  (a+ib)(x+iy) = ax-by + i(xb+ay)
        if( (Math.abs(c.real) < Math.eps) && (Math.abs(c.imaginary) < Math.eps) ){
            this.real = Infinity;
            this.imaginary = Infinity;

            return;
        }

        denom = c.real*c.real + c.imaginary*c.imaginary;

        this.real = (this.real*c.real + this.imaginary*c.imaginary)/denom;
        this.imaginary = (this.imaginary*c.real - this.real*c.imaginary)/denom;
    }
};


/**
 * @description
 * JXG.C is the complex number (name)space. It provides functions to calculate with
 * complex numbers (defined in {@link JXG.Complex}). With this namespace you don't have to modify
 * your existing complex numbers, e.g. to add two complex numbers:
 * <pre class="code">   var z1 = new JXG.Complex(1, 0);
 *    var z2 = new JXG.Complex(0, 1);
 *    z = JXG.C.add(z1, z1);</pre>
 * z1 and z2 here remain unmodified. With the object oriented approach above this
 * section the code would look like:
 * <pre class="code">   var z1 = new JXG.Complex(1, 0);
 *    var z2 = new JXG.Complex(0, 1);
 *    var z = new JXG.Complex(z1);
 *    z.add(z2);</pre>
 * @namespace Namespace for the complex number arithmetic functions.
 */
JXG.C = {};

/**
 * Add two (complex) numbers z1 and z2 and return the result as a (complex) number.
 * @param z1 Summand
 * @param z2 Summand
 * @return A complex number equal to the sum of the given parameters.
 */
JXG.C.add = function(/** JXG.Complex,number */ z1, /** JXG.Complex,number */ z2) /** JXG.Complex */{
    var z = new JXG.Complex(z1);
    z.add(z2);
    return z;
};

/**
 * Subtract two (complex) numbers z1 and z2 and return the result as a (complex) number.
 * @param z1 Minuend
 * @param z2 Subtrahend
 * @return A complex number equal to the difference of the given parameters.
 */
JXG.C.sub = function(/** JXG.Complex,number */ z1, /** JXG.Complex,number */ z2) /** JXG.Complex */{
    var z = new JXG.Complex(z1);
    z.sub(z2);
    return z;
};

/**
 * Multiply two (complex) numbers z1 and z2 and return the result as a (complex) number.
 * @param z1 Factor
 * @param z2 Factor
 * @return A complex number equal to the product of the given parameters.
 */
JXG.C.mult = function(/** JXG.Complex,number */ z1, /** JXG.Complex,number */ z2) /** JXG.Complex */{
    var z = new JXG.Complex(z1);
    z.mult(z2);
    return z;
};

/**
 * Divide two (complex) numbers z1 and z2 and return the result as a (complex) number.
 * @param z1 Dividend
 * @param z2 Divisor
 * @return A complex number equal to the quotient of the given parameters.
 */
JXG.C.div = function(/** JXG.Complex,number */ z1, /** JXG.Complex,number */ z2) /** JXG.Complex */{
    var z = new JXG.Complex(z1);
    z.div(z2);
    return z;
};

/* 
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.

*/

/** 
 * @fileoverview AbstractRenderer is the base class for all renderers. This class is subject
 * to change the next time that only a few methods are strictly required for a special
 * renderer derived from this class to have a working special renderer. Of course other
 * methods may be overwritten for performance reasons, too.
 */

/**
 * Constructs a new AbstractRenderer object.
 * @class The AbstractRenderer is a base class that should be considered what in other languages
 * is called an abstract class, even though no such thing really exists in JavaScript. All members
 * essential for a renderer are defined here.
 * @constructor
 * @see JXG.SVGRenderer
 * @see JXG.VMLRenderer
 * @private
 */
JXG.AbstractRenderer = function() {
	/**
	 * TODO Needs description
	 * @type number
	 * @default 8
	 * @private
	 */
    this.vOffsetText = 8;
    
    /**
     * If true the visual properties of the elements are updated on every update.
     * @type boolean
     * @default true 
     */
    this.enhancedRendering = true;
};


/* ************************** 
 *    Point related stuff
 * **************************/

/**
 * Draws a point on the canvas.
 * @param el Reference to a point object, that has to be drawn.
 * @see JXG.Point
 * @see #updatePoint
 */
JXG.AbstractRenderer.prototype.drawPoint = function(/** JXG.Point */ el) {
    var node,
        f = el.visProp['face'];
        
    if(f == 'cross' || f == 'x') { // x
        node = this.createPrimitive('path',el.id);
        this.appendChildPrimitive(node,el.layer);
        this.appendNodesToElement(el, 'path');
    }
    else if(f == 'circle' || f == 'o') { // circle
        node = this.createPrimitive('circle',el.id);
        this.appendChildPrimitive(node,el.layer);
        this.appendNodesToElement(el, 'circle');
    }
    else if(f == 'square' || f == '[]') { // rectangle
        node = this.createPrimitive('rect',el.id);
        this.appendChildPrimitive(node,el.layer);
        this.appendNodesToElement(el, 'rect');
    }
    else if(f == 'plus' || f == '+') { // +
        node = this.createPrimitive('path',el.id);
        this.appendChildPrimitive(node,el.layer);  
        this.appendNodesToElement(el, 'path');
    }
    else if(f == 'diamond' || f == '<>') {
        node = this.createPrimitive('path',el.id);
        this.appendChildPrimitive(node,el.layer);  
        this.appendNodesToElement(el, 'path');
    }
    else if(f == 'triangleup' || f == 'a') {
        node = this.createPrimitive('path',el.id);
        this.appendChildPrimitive(node,el.layer);  
        this.appendNodesToElement(el, 'path');    
    }
    else if(f == 'triangledown' || f == 'v') {
        node = this.createPrimitive('path',el.id);
        this.appendChildPrimitive(node,el.layer);  
        this.appendNodesToElement(el, 'path');
    }        
    else if(f == 'triangleleft' || f == '<') {
        node = this.createPrimitive('path',el.id);
        this.appendChildPrimitive(node,el.layer);  
        this.appendNodesToElement(el, 'path');    
    }
    else if(f == 'triangleright' || f == '>') {
        node = this.createPrimitive('path',el.id);
        this.appendChildPrimitive(node,el.layer);  
        this.appendNodesToElement(el, 'path');    
    }
    el.rendNode = node;
    
    this.setObjectStrokeWidth(el,el.visProp['strokeWidth']);
    this.setObjectStrokeColor(el,el.visProp['strokeColor'],el.visProp['strokeOpacity']);
    this.setObjectFillColor(el,el.visProp['fillColor'],el.visProp['fillOpacity']);	

    this.updatePoint(el);
};
   
/**
 * Updates visual appearance of the renderer element assigned to the given line.
 * @param el Reference to a point object, that has to be updated.
 * @see JXG.Point
 * @see #drawPoint
 * @see #changePointStyle
 */
JXG.AbstractRenderer.prototype.updatePoint = function(/** JXG.Point */ el) {
    var size = el.visProp['size'],
    	f = el.visProp['face'];
    if (isNaN(el.coords.scrCoords[2]) || isNaN(el.coords.scrCoords[1])) return;
    
    if (this.enhancedRendering) {
        if (!el.visProp['draft']) {
            this.setObjectStrokeWidth(el,el.visProp['strokeWidth']);
            this.setObjectStrokeColor(el,el.visProp['strokeColor'],el.visProp['strokeOpacity']);
            this.setObjectFillColor(el,el.visProp['fillColor'],el.visProp['fillOpacity']);
        } else {
            this.setDraft(el);
        }
    }
    // Zoom does not work for traces.
    size *= ((!el.board || !el.board.options.point.zoom)?1.0:Math.sqrt(el.board.zoomX*el.board.zoomY));
    
    if(f == 'cross' || f == 'x') { // x
        this.updatePathPrimitive(el.rendNode, this.updatePathStringPoint(el, size,'x'), el.board); 
    }
    else if(f == 'circle' || f == 'o') { // circle
        this.updateCirclePrimitive(el.rendNode,el.coords.scrCoords[1], el.coords.scrCoords[2], size+1);            
    }
    else if(f == 'square' || f == '[]') { // rectangle
        this.updateRectPrimitive(el.rendNode,
                el.coords.scrCoords[1]-size, el.coords.scrCoords[2]-size, size*2, size*2);
    }
    else if(f == 'plus' || f == '+') { // +
        this.updatePathPrimitive(el.rendNode, this.updatePathStringPoint(el,size,'+'), el.board); 
    }
    else if(f == 'diamond' || f == '<>') { // diamond
        this.updatePathPrimitive(el.rendNode, this.updatePathStringPoint(el,size,'diamond'), el.board); 
    }  
    else if(f == 'triangleup' || f == 'a') { // triangleUp
        this.updatePathPrimitive(el.rendNode, this.updatePathStringPoint(el,size,'A'), el.board); 
    } 
    else if(f == 'triangledown' || f == 'v') { // triangleDown
        this.updatePathPrimitive(el.rendNode, this.updatePathStringPoint(el,size,'v'), el.board); 
    } 
    else if(f == 'triangleleft' || f == '<') { // triangleLeft
        this.updatePathPrimitive(el.rendNode, this.updatePathStringPoint(el,size,'<'), el.board); 
    }  
    else if(f == 'triangleright' || f == '>') { // triangleRight
        this.updatePathPrimitive(el.rendNode, this.updatePathStringPoint(el,size,'>'), el.board); 
    }        
    this.setShadow(el);
};

/**
 * Changes the style of a point that already exists on the canvas. This is required because
 * the point styles differ in what elements have to be drawn, e.g. if the point is marked by
 * a x or a + two lines are drawn, if it's marked by spot a circle is drawn. This method removes
 * the old renderer element(s) and creates the new one(s).
 * @param el Reference to a point object, that has to be updated.
 * @see JXG.Point
 * @see #updatePoint
 */
JXG.AbstractRenderer.prototype.changePointStyle = function(/** JXG.Point */el) {
    var node = this.getElementById(el.id);
    if(node != null) {
        this.remove(node);
    }
    this.drawPoint(el);
    JXG.clearVisPropOld(el);

    if(!el.visProp['visible']) {
        this.hide(el);
    }
    if(el.visProp['draft']) {
        this.setDraft(el);
    }
};


/* ************************** 
 *    Line related stuff
 * **************************/


/**
 * Draws a line on the canvas.
 * @param {JXG.Line} el Reference to a line object, that has to be drawn.
 * @see JXG.Line
 * @see #updateLine
 * @see #calcStraight
 */
JXG.AbstractRenderer.prototype.drawLine = function(el) { 
    var node = this.createPrimitive('line',el.id);
    this.appendChildPrimitive(node,el.layer);
    this.appendNodesToElement(el,'lines');

    this.updateLine(el);
};

/**
 * Updates visual appearance of the renderer element assigned to the given line.
 * @param el Reference to the line object that has to be updated.
 * @see JXG.Line
 * @see #drawLine
 * @see #calcStraight
 */
JXG.AbstractRenderer.prototype.updateLine = function(/** JXG.Line */ el) {
    var screenCoords1 = new JXG.Coords(JXG.COORDS_BY_USER, el.point1.coords.usrCoords, el.board),
        screenCoords2 = new JXG.Coords(JXG.COORDS_BY_USER, el.point2.coords.usrCoords, el.board),
        ax, ay, bx, by, beta, sgn, x, y, m;
        
    //if(el.visProp['straightFirst'] || el.visProp['straightLast']) {
       this.calcStraight(el,screenCoords1,screenCoords2); 
    //} 
    this.updateLinePrimitive(el.rendNode,screenCoords1.scrCoords[1],screenCoords1.scrCoords[2],
            screenCoords2.scrCoords[1],screenCoords2.scrCoords[2],el.board);

    // Update the image which is connected to the line:
    if (el.image!=null) {
        ax = screenCoords1.scrCoords[1];
        ay = screenCoords1.scrCoords[2];
        bx = screenCoords2.scrCoords[1];
        by = screenCoords2.scrCoords[2];
        //beta;                                              // ???
        sgn = (bx-ax>=0)?1:-1;
        if (Math.abs(bx-ax)>0.0000001) {
            beta = Math.atan2(by-ay,bx-ax)+ ((sgn<0)?Math.PI:0);  
        } else {
            beta = ((by-ay>0)?0.5:-0.5)*Math.PI;
        }
        x = 250; //ax;
        y = 256; //ay;//+el.image.size[1]*0.5;
        m = [
                 [1,                                    0,             0],
                 [x*(1-Math.cos(beta))+y*Math.sin(beta),Math.cos(beta),-Math.sin(beta)],
                 [y*(1-Math.cos(beta))-x*Math.sin(beta),Math.sin(beta), Math.cos(beta)]
                ];
        el.imageTransformMatrix = m;
    }
    this.makeArrows(el);
    
    if (this.enhancedRendering) {
        if (!el.visProp['draft']) {
            this.setObjectStrokeWidth(el,el.visProp['strokeWidth']);
            this.setObjectStrokeColor(el,el.visProp['strokeColor'],el.visProp['strokeOpacity']);
            this.setDashStyle(el,el.visProp);    
            this.setShadow(el);
        } else {
            this.setDraft(el);
        }
    }     
};

/**
 * Calculates drawing start and end point for a line. A segment is only drawn from start to end point, a straight line
 * is drawn until it meets the boards boundaries.
 * @param el Reference to a line object, that needs calculation of start and end point.
 * @param point1 Coordinates of the point where line drawing begins.
 * @param point2 Coordinates of the point where line drawing ends.
 * @see JXG.Line
 * @see #drawLine
 * @see #updateLine
 */
JXG.AbstractRenderer.prototype.calcStraight = function(/** JXG.Line */ el, /** JXG.Coords */ point1, /** JXG.Coords */ point2) {
    var takePoint1, takePoint2, intersect1, intersect2, straightFirst, straightLast, 
        b, c, s, i, j, p1, p2;
    
    //b = el.board.algebra;
    straightFirst = el.visProp['straightFirst'];
    straightLast  = el.visProp['straightLast'];

/*
    if (Math.abs(point1.scrCoords[0])<b.eps||Math.abs(point2.scrCoords[0])<b.eps) {
        straightFirst = true;
        straightLast  = true;
    }
*/    
    // If one of the point is an ideal point in homogeneous coordinates
    // drawing of line segments or rays are not possible. 
    if (Math.abs(point1.scrCoords[0])<JXG.Math.eps) {
        straightFirst = true;
    }
    if (Math.abs(point2.scrCoords[0])<JXG.Math.eps) {
        straightLast  = true;
    }

    if ( !straightFirst && !straightLast ) {  // Do nothing in case of line segments (inside or outside of the board)
        return;
    }
    
    // Compute the stdform of the line in screen coordinates.
    c = [];
    c[0] = el.stdform[0] - 
           el.stdform[1]*el.board.origin.scrCoords[1]/el.board.stretchX+
           el.stdform[2]*el.board.origin.scrCoords[2]/el.board.stretchY;
    c[1] = el.stdform[1]/el.board.stretchX;
    c[2] = el.stdform[2]/(-el.board.stretchY);

    if (isNaN(c[0]+c[1]+c[2])) return; // p1=p2
    
    // Intersect the line with the four borders of the board.
    s = [];
    s[0] = JXG.Math.crossProduct(c,[0,0,1]);  // top
    s[1] = JXG.Math.crossProduct(c,[0,1,0]);  // left
    s[2] = JXG.Math.crossProduct(c,[-el.board.canvasHeight,0,1]);  // bottom
    s[3] = JXG.Math.crossProduct(c,[-el.board.canvasWidth,1,0]);   // right

    // Normalize the intersections 
    for (i=0;i<4;i++) {
        if (Math.abs(s[i][0])>JXG.Math.eps) {
            for (j=2;j>0;j--) {
                s[i][j] /= s[i][0];
            }
            s[i][0] = 1.0;
        }
    }
    
    takePoint1 = false;
    takePoint2 = false;
    if (!straightFirst &&    // Line starts at point1 and point2 is inside the board
            point1.scrCoords[1]>=0.0 && point1.scrCoords[1]<=el.board.canvasWidth &&
            point1.scrCoords[2]>=0.0 && point1.scrCoords[2]<=el.board.canvasHeight) {
        takePoint1 = true;
    }
    if (!straightLast &&    // Line ends at point2 and point2 is inside the board
            point2.scrCoords[1]>=0.0 && point2.scrCoords[1]<=el.board.canvasWidth &&
            point2.scrCoords[2]>=0.0 && point2.scrCoords[2]<=el.board.canvasHeight) {
        takePoint2 = true;
    }

    if (Math.abs(s[1][0])<JXG.Math.eps) {                  // line is parallel to "left", take "top" and "bottom"
        intersect1 = s[0];                          // top
        intersect2 = s[2];                          // bottom
    } else if (Math.abs(s[0][0])<JXG.Math.eps) {           // line is parallel to "top", take "left" and "right"
        intersect1 = s[1];                          // left
        intersect2 = s[3];                          // right
    } else if (s[1][2]<0) {                         // left intersection out of board (above)
        intersect1 = s[0];                          // top
        if (s[3][2]>el.board.canvasHeight) {        // right intersection out of board (below)
            intersect2 = s[2];                      // bottom
        } else {
            intersect2 = s[3];                      // right
        }
    } else if (s[1][2]>el.board.canvasHeight) {     // left intersection out of board (below)
        intersect1 = s[2];                          // bottom
        if (s[3][2]<0) {                            // right intersection out of board (above)
            intersect2 = s[0];                      // top
        } else {
            intersect2 = s[3];                      // right
        }
    } else {
        intersect1 = s[1];                          // left
        if (s[3][2]<0) {                            // right intersection out of board (above)
            intersect2 = s[0];                      // top
        } else if (s[3][2]>el.board.canvasHeight) { // right intersection out of board (below)
            intersect2 = s[2];                      // bottom
        } else {
            intersect2 = s[3];                      // right
        }
    }
    
    intersect1 = new JXG.Coords(JXG.COORDS_BY_SCREEN, intersect1.slice(1), el.board);
    intersect2 = new JXG.Coords(JXG.COORDS_BY_SCREEN, intersect2.slice(1), el.board);
 
    if (!takePoint1) {
        if (!takePoint2) {                // Two border intersection points are used
            if (this.isSameDirection(point1, point2, intersect1)) {
                if (!this.isSameDirection(point1, point2, intersect2)) {
                    p2 = intersect1;
                    p1 = intersect2;
                } else {
                    if (el.board.algebra.affineDistance(point2.usrCoords,intersect1.usrCoords)<el.board.algebra.affineDistance(point2.usrCoords,intersect2.usrCoords)) {
                        p1 = intersect1;
                        p2 = intersect2;
                    } else {
                        p2 = intersect1;
                        p1 = intersect2;
                    }
                }
            } else {
                if (this.isSameDirection(point1, point2, intersect2)) {
                    p1 = intersect1;
                    p2 = intersect2;
                } else {
                    if (el.board.algebra.affineDistance(point2.usrCoords,intersect1.usrCoords)<el.board.algebra.affineDistance(point2.usrCoords,intersect2.usrCoords)) {
                        p2 = intersect1;
                        p1 = intersect2;
                    } else {
                        p1 = intersect1;
                        p2 = intersect2;
                    }
                }
            }
        } else {                          // Instead of point1 the border intersection is taken
            if (this.isSameDirection(point2, point1, intersect1)) {
                p1 = intersect1;
            } else {
                p1 = intersect2;
            }
        }
    } else {
        if (!takePoint2) {                // Instead of point2 the border intersection is taken
            if (this.isSameDirection(point1, point2, intersect1)) {
                p2 = intersect1;
            } else {
                p2 = intersect2;
            }
        }
    }

    if (p1) point1.setCoordinates(JXG.COORDS_BY_USER, p1.usrCoords.slice(1));
    if (p2) point2.setCoordinates(JXG.COORDS_BY_USER, p2.usrCoords.slice(1));
};

/**
* If you're looking from point "start" towards point "s" and can see the point "p", true is returned. Otherwise false.
* @param start The point you're standing on.
* @param p The point in which direction you're looking.
* @param s The point that should be visible.
* @return True, if from start the point p is in the same direction as s is, that means s-start = k*(p-start) with k>=0.
* @private
*/
JXG.AbstractRenderer.prototype.isSameDirection = function(/** JXG.Coords */ start, /** JXG.Coords */ p, /** JXG.Coords */ s) /** boolean */ {
    var dx, dy, sx, sy;
    
    dx = p.usrCoords[1]-start.usrCoords[1];
    dy = p.usrCoords[2]-start.usrCoords[2];

    sx = s.usrCoords[1]-start.usrCoords[1];
    sy = s.usrCoords[2]-start.usrCoords[2];

    if (Math.abs(dx)<JXG.Math.eps) dx=0;
    if (Math.abs(dy)<JXG.Math.eps) dy=0;
    if (Math.abs(sx)<JXG.Math.eps) sx=0;
    if (Math.abs(sy)<JXG.Math.eps) sy=0;

    if (dx>=0&&sx>=0) {
        if ((dy>=0&&sy>=0) || (dy<=0&&sy<=0)) { return true; }
    } else if (dx<=0&&sx<=0){
        if ((dy>=0&&sy>=0) || (dy<=0&&sy<=0)) { return true; }        
    }

    return false;
};

/**
 * Update ticks on a line. This method is only a stub and is by now implemented only in the special renderers.
 * @param axis Reference of an line object, thats ticks have to be updated.
 * @param dxMaj Number of pixels a major tick counts in x direction.
 * @param dyMaj Number of pixels a major tick counts in y direction.
 * @param dxMin Number of pixels a minor tick counts in x direction.
 * @param dyMin Number of pixels a minor tick counts in y direction.
 * @see JXG.Line
 */
JXG.AbstractRenderer.prototype.updateTicks = function(/** JXG.Line */ axis, /** number */ dxMaj,
                                                      /** number */ dyMaj, /** number */ dxMin, /** number */ dyMin)
{ };

/**
 * Removes all ticks from an axis.
 * @param axis Reference of an line object, that's ticks have to be removed.
 * @deprecated
 * @see JXG.Line
 * @see #upateTicks
 */
JXG.AbstractRenderer.prototype.removeTicks = function(/** JXG.Line */ axis) {
    var ticks = this.getElementById(axis.id+'_ticks');
    this.remove(ticks);
};


/**
 * Creates a rendering node for an arrow on the board.
 * @param el Reference to an arrow object, that has to be drawn.
 * @see Arrow
 * @see JXG.Line
 * @see #updateArrow
 */
JXG.AbstractRenderer.prototype.drawArrow = function(/** JXG.Line */ el) {
    var node = this.createPrimitive('line',el.id);
    this.setObjectStrokeWidth(el,el.visProp['strokeWidth']); // ?
    this.setObjectStrokeColor(el,el.visProp['strokeColor'],el.visProp['strokeOpacity']); // ?
    this.setObjectFillColor(el,el.visProp['fillColor'],el.visProp['fillOpacity']); // ?
    this.setDashStyle(el,el.visProp); // ?
    this.makeArrow(node,el);
    this.appendChildPrimitive(node,el.layer);
    this.appendNodesToElement(el,'lines');

    this.updateArrow(el);
};

/**
 * Updates properties of an arrow that already exists on the canvas.
 * @param el Reference to an arrow object, that has to be updated.
 * @see Arrow
 * @see JXG.Line
 * @see #drawArrow
 */
JXG.AbstractRenderer.prototype.updateArrow = function(/** JXG.Line */ el) {
    if (this.enhancedRendering) {
        if (!el.visProp['draft']) {
            this.setObjectStrokeWidth(el,el.visProp['strokeWidth']);
            this.setObjectStrokeColor(el,el.visProp['strokeColor'],el.visProp['strokeOpacity']);
            this.setObjectFillColor(el,el.visProp['fillColor'],el.visProp['fillOpacity']);
            this.setShadow(el);
            this.setDashStyle(el,el.visProp);
        } else {
            this.setDraft(el);
        }
    }
    this.updateLinePrimitive(el.rendNode,el.point1.coords.scrCoords[1],el.point1.coords.scrCoords[2],
        el.point2.coords.scrCoords[1],el.point2.coords.scrCoords[2],el.board);
};


/* ************************** 
 *    Curve related stuff
 * **************************/

/**
 * Draws a graph on the canvas.
 * @param {JXG.Curve} el Reference to a graph object, that has to be plotted.
 * @see JXG.Curve
 * @see #updateCurve
 */
JXG.AbstractRenderer.prototype.drawCurve = function(el) { 
    var node = this.createPrimitive('path',el.id);
    
    //node.setAttributeNS(null, 'stroke-linejoin', 'round');
    this.appendChildPrimitive(node,el.layer);
    this.appendNodesToElement(el,'path');
    this.setObjectStrokeWidth(el,el.visProp['strokeWidth']); // ?
    this.setObjectStrokeColor(el,el.visProp['strokeColor'],el.visProp['strokeOpacity']); // ?
    this.setObjectFillColor(el,el.visProp['fillColor'],el.visProp['fillOpacity']); // ?
    this.setDashStyle(el,el.visProp); // ?
    this.updateCurve(el);
};

/**
 * Updates visual appearance of the renderer element assigned to the given curve.
 * @param el Reference to a curve object, that has to be updated.
 * @see JXG.Curve
 * @see #drawCurve
 */
JXG.AbstractRenderer.prototype.updateCurve = function(/** JXG.Curve */ el) {
    if (this.enhancedRendering) {
        if (!el.visProp['draft']) {
            this.setObjectStrokeWidth(el,el.visProp['strokeWidth']);
            this.setObjectStrokeColor(el,el.visProp['strokeColor'],el.visProp['strokeOpacity']);
            this.setObjectFillColor(el,el.visProp['fillColor'],el.visProp['fillOpacity']);
            this.setDashStyle(el,el.visProp);
            this.setShadow(el);
        } else {
            this.setDraft(el);
        }
    }
    this.updatePathPrimitive(el.rendNode,this.updatePathStringPrimitive(el),el.board);
};


/* ************************** 
 *    Circle related stuff
 * **************************/

/**
 * Draws a circle on the canvas.
 * @param el Reference to a circle object, that has to be drawn.
 * @see JXG.Circle
 * @see #updateCircle
 */
JXG.AbstractRenderer.prototype.drawCircle = function(/** JXG.Circle */ el) { 
    var node = this.createPrimitive('ellipse',el.id);
    this.appendChildPrimitive(node,el.layer);
    this.appendNodesToElement(el,'ellipse'); 
    
    this.updateCircle(el);
};

/**
 * Updates visual appearance of a given circle on the board.
 * @param {JXG.Circle} el Reference to a circle object, that has to be updated.
 * @see JXG.Circle
 * @see #drawCircle
 */
JXG.AbstractRenderer.prototype.updateCircle = function(el) {
    if (this.enhancedRendering) {
        if (!el.visProp['draft']) {
            this.setObjectStrokeWidth(el,el.visProp['strokeWidth']);
            this.setObjectStrokeColor(el,el.visProp['strokeColor'],el.visProp['strokeOpacity']);
            this.setObjectFillColor(el,el.visProp['fillColor'],el.visProp['fillOpacity']);
            this.setDashStyle(el,el.visProp);
            this.setShadow(el);
        } else {
            this.setDraft(el);
        }
    }
    // Radius umrechnen:
    var radius = el.Radius();
    if (radius>0.0 && !isNaN(el.midpoint.coords.scrCoords[1]+el.midpoint.coords.scrCoords[2]) ) {
        this.updateEllipsePrimitive(el.rendNode,el.midpoint.coords.scrCoords[1],el.midpoint.coords.scrCoords[2],
            (radius * el.board.stretchX),(radius * el.board.stretchY));
    }
};
    

/* ************************** 
 *   Polygon related stuff
 * **************************/

/**
 * Draws a polygon on the canvas.
 * @param el Reference to a Polygon object, that is to be drawn.
 * @see JXG.Polygon
 * @see #updatePolygon
 */
JXG.AbstractRenderer.prototype.drawPolygon = function(/** JXG.Polygon */ el) { 
    var node = this.createPrimitive('polygon',el.id);
    el.visProp['fillOpacity'] = 0.3;
    //el.visProp['strokeColor'] = 'none';
    //this.setObjectFillColor(el,el.visProp['fillColor'],el.visProp['fillOpacity']);
    this.appendChildPrimitive(node,el.layer);
    this.appendNodesToElement(el,'polygon');
    this.updatePolygon(el);
};
    
/**
 * Updates properties of a polygon's rendering node.
 * @param {} el Reference to a polygon object, that has to be updated.
 * @see JXG.Polygon
 * @see #drawPolygon
 */
JXG.AbstractRenderer.prototype.updatePolygon = function(/** JXG.Polygon */ el) { 
    if (this.enhancedRendering) {
        if (!el.visProp['draft']) {
            this.setObjectStrokeWidth(el,el.visProp['strokeWidth']);
            this.setObjectFillColor(el,el.visProp['fillColor'],el.visProp['fillOpacity']);
            this.setShadow(el); 
        } else {
            this.setDraft(el);
        }
    }

    this.updatePolygonePrimitive(el.rendNode,el);
};


/* ************************** 
 *    Arc related stuff
 * **************************/

/**
 * Draws an arc on the canvas; This method is a stub and has to be implemented by the special renderers.
 * @param arc Reference to an arc object, that has to be drawn.
 * @see JXG.Arc
 * @see #updateArc
 */
JXG.AbstractRenderer.prototype.drawArc = function(/** JXG.Arc */ arc) { };

/**
 * Updates properties of an arc; This method is a stub and has to be implemented by the special renderers.
 * @param arc Reference to an arc object, that has to be updated.
 * @see JXG.Arc
 * @see #drawArc
 */
JXG.AbstractRenderer.prototype.updateArc = function(/** JXG.Arc */ el) { };


/* ************************** 
 *    Text related stuff
 * **************************/

/**
 * Puts a text node onto the board.
 * @param text Reference to an text object, that has to be drawn
 * @see JXG.Text
 * @see #updateText
 * @see #updateTextStyle
 */
JXG.AbstractRenderer.prototype.drawText = function(/** JXG.Text */ el) { 
    var node;
    if (el.display=='html') {
        node = this.container.ownerDocument.createElement('div');
        node.style.position = 'absolute';
        node.style.fontSize = el.board.fontSize + 'px';  
        node.style.color = el.visProp['strokeColor'];
        node.className = 'JXGtext';
        node.style.zIndex = '10';      
        this.container.appendChild(node);
        node.setAttribute('id', el.id);
    } else {
        node = this.drawInternalText(el);
    }
    el.rendNode = node;
    el.htmlStr = '';
    this.updateText(el);
};

JXG.AbstractRenderer.prototype.drawInternalText = function(el) {};


/**
 * Updates GeometryElement properties of an already existing text element.
 * @param el Reference to an text object, that has to be updated.
 * @see JXG.Text
 * @see #drawText
 * @see #updateTextStyle
 */
JXG.AbstractRenderer.prototype.updateText = function(/** JXG.Text */ el) { 
    // Update onky objects that are visible.
    if (el.visProp['visible'] == false) return;
    if (isNaN(el.coords.scrCoords[1]+el.coords.scrCoords[2])) return;
    this.updateTextStyle(el);
    if (el.display=='html') {
        el.rendNode.style.left = (el.coords.scrCoords[1])+'px'; 
        el.rendNode.style.top = (el.coords.scrCoords[2] - this.vOffsetText)+'px'; 
        el.updateText();
        if (el.htmlStr!= el.plaintextStr) {
            el.rendNode.innerHTML = el.plaintextStr;
            if (el.board.options.text.useASCIIMathML) {
                AMprocessNode(el.rendNode,false);
            }
            el.htmlStr = el.plaintextStr;
        }
    } else {
        this.updateInternalText(el);
    }
};

JXG.AbstractRenderer.prototype.updateInternalText = function(el) {};

/**
 * Updates CSS style properties of a text node.
 * @param el Reference to the text object, that has to be updated.
 * @see JXG.Text
 * @see #drawText
 * @see #updateText
 */
JXG.AbstractRenderer.prototype.updateTextStyle = function(/** JXG.Text */ el) { 
    var fs;
    if (el.visProp['fontSize']) {
        if (typeof el.visProp['fontSize'] == 'function') {
            fs = el.visProp['fontSize']();
            el.rendNode.style.fontSize = (fs>0?fs:0); 
        } else {
            el.rendNode.style.fontSize = (el.visProp['fontSize']); 
        }
    }
};


/* ************************** 
 *    Angle related stuff
 * **************************/

/**
 * Draft method for special renderers to draw an angle.
 * @param angle Reference to an angle object, that has to be drawn.
 * @see JXG.Angle
 * @see #updateAngle
 */
JXG.AbstractRenderer.prototype.drawAngle = function(/** JXG.Angle */ angle) { };

/**
 * Update method draft for updating the properties of an angle.
 * @param angle Reference to an angle object.
 * @see JXG.Angle
 * @see #drawAngle
 */
JXG.AbstractRenderer.prototype.updateAngle = function(/** JXG.Angle */ angle) { };


/* ************************** 
 *    Image related stuff
 * **************************/

/**
 * Draws an image on the canvas; This is just a template, has to be implemented by special renderers.
 * @param image Reference to an image object, that has to be drawn.
 * @see JXG.Image
 * @see #updateImage
 */
JXG.AbstractRenderer.prototype.drawImage = function(/** JXG.Image */ image) { };

/**
 * Updates the properties of an Image element.
 * @param el Reference to an image object, that has to be updated.
 * @see JXG.Image
 * @see #drawImage
 */
JXG.AbstractRenderer.prototype.updateImage = function(/** JXG.Image */ el) { 
    this.updateRectPrimitive(el.rendNode,el.coords.scrCoords[1],el.coords.scrCoords[2]-el.size[1],
        el.size[0],el.size[1]);
        
    if (el.parent != null) {
        this.transformImageParent(el,el.parent.imageTransformMatrix);
    } else {
        this.transformImageParent(el); // Transforms are cleared
    }
    this.transformImage(el,el.transformations);    
};


/* ************************** 
 *    Grid stuff
 * **************************/

/**
 * Creates a grid on the board, i.e. light helper lines to support the user on creating and manipulating a construction.
 * @param board Board on which the grid is drawn.
 * @see #removeGrid
 */
JXG.AbstractRenderer.prototype.drawGrid = function(/** JXG.Board */ board) { 
    var gridX = board.gridX,
        gridY = board.gridY,
        k = new JXG.Coords(JXG.COORDS_BY_SCREEN, [0,0], board),
        k2 = new JXG.Coords(JXG.COORDS_BY_SCREEN, [board.canvasWidth, board.canvasHeight], board),
        tmp = Math.ceil(k.usrCoords[1]),
        j = 0,
        i, j2, l, l2,
        gx, gy, topLeft, bottomRight, node2,
        el, eltmp;
        
    board.hasGrid = true;

    for(i = 0; i <= gridX+1; i++) {
        if(tmp-i/gridX < k.usrCoords[1]) {
            j = i-1;
            break;
        }
    }

    tmp = Math.floor(k2.usrCoords[1]);
    j2 = 0;
    for(i = 0; i <= gridX+1; i++) {
        if(tmp+i/gridX > k2.usrCoords[1]) {
            j2 = i-1;
            break;
        }
    } 

    tmp = Math.ceil(k2.usrCoords[2]);
    l2 = 0;
    for(i = 0; i <= gridY+1; i++) {
        if(tmp-i/gridY < k2.usrCoords[2]) {
            l2 = i-1;
            break;
        }
    }

    tmp = Math.floor(k.usrCoords[2]);
    l = 0;
    for(i = 0; i <= gridY+1; i++) {
        if(tmp+i/gridY > k.usrCoords[2]) {
            l = i-1;
            break;
        }
    }

    gx = Math.round((1.0/gridX)*board.stretchX);
    gy = Math.round((1.0/gridY)*board.stretchY);

    topLeft = new JXG.Coords(JXG.COORDS_BY_USER, 
                             [Math.ceil(k.usrCoords[1])-j/gridX, Math.floor(k.usrCoords[2])+l/gridY],
                             board);
    bottomRight = new JXG.Coords(JXG.COORDS_BY_USER,
                                 [Math.floor(k2.usrCoords[1])+j2/gridX, Math.ceil(k2.usrCoords[2])-l2/gridY],
                                 board);
                                     
    node2 = this.drawVerticalGrid(topLeft, bottomRight, gx, board);
    this.appendChildPrimitive(node2, board.options.layer['grid']);
    if(!board.snapToGrid) {
        el = new Object();
        el.rendNode = node2;
        el.elementClass = JXG.OBJECT_CLASS_LINE;
        el.id = "gridx";
        JXG.clearVisPropOld(el);
        this.setObjectStrokeColor(el, board.gridColor, board.gridOpacity);
    }
    else {
        el = new Object();
        el.rendNode = node2;
        el.elementClass = JXG.OBJECT_CLASS_LINE;
        el.id = "gridx";        
        JXG.clearVisPropOld(el);
        this.setObjectStrokeColor(el, '#FF8080', 0.5); //board.gridOpacity);    
    }
    this.setPropertyPrimitive(node2,'stroke-width', '0.4px');  
    if(board.gridDash) {
        this.setGridDash("gridx"); 
    }

    node2 = this.drawHorizontalGrid(topLeft, bottomRight, gy, board);
    this.appendChildPrimitive(node2, board.options.layer['grid']); // Attention layer=1
    if(!board.snapToGrid) {
        el = new Object();
        el.rendNode = node2;
        el.elementClass = JXG.OBJECT_CLASS_LINE;
        el.id = "gridy";   
        JXG.clearVisPropOld(el);
        this.setObjectStrokeColor(el, board.gridColor, board.gridOpacity);
    }
    else {
        el = new Object();
        el.rendNode = node2;
        el.elementClass = JXG.OBJECT_CLASS_LINE;
        el.id = "gridy";        
        JXG.clearVisPropOld(el);
        this.setObjectStrokeColor(el, '#FF8080', 0.5); //board.gridOpacity);    
    }
    this.setPropertyPrimitive(node2,'stroke-width', '0.4px');  
    if(board.gridDash) {
        this.setGridDash("gridy"); 
    }
   
};

/**
 * Remove the grid from the given board; This is a template that has to be implemented by the renderer itself.
 * @param board Board from which the grid is removed.
 * @see #drawGrid
 */
JXG.AbstractRenderer.prototype.removeGrid = function(/** JXG.Board */ board) {
    var c = document.getElementById('gridx');
    this.remove(c);

    c = document.getElementById('gridy');
    this.remove(c);

    board.hasGrid = false;
/*
    var c = this.layer[board.options.layer['grid']];
    board.hasGrid = false;
    while (c.childNodes.length>0) {
        c.removeChild(c.firstChild);
    }
 */
};
 


/* ************************** 
 *  general element helpers
 * **************************/

/**
 * Hides an element on the canvas; Only a stub, requires implementation in the derived renderer.
 * @param obj Reference to the geometry element that has to disappear.
 * @see #show
 */
JXG.AbstractRenderer.prototype.hide = function(/** JXG.GeometryElement */ obj) { };

/**
 * Shows a hidden element on the canvas; Only a stub, requires implementation in the derived renderer.
 * @param obj Reference to the object that has to appear.
 * @see #hide
 */
JXG.AbstractRenderer.prototype.show = function(/** JXG.GeometryElement */ obj) { };

/**
 * Sets an element's stroke width.
 * @param el Reference to the geometry element.
 * @param width The new stroke width to be assigned to the element.
 */
JXG.AbstractRenderer.prototype.setObjectStrokeWidth = function(/** JXG.GeometryElement */ el, /** number */ width) { };

/**
 * Changes an objects stroke color to the given color.
 * @param obj Reference of the {@link JXG.GeometryElement} that gets a new stroke color.
 * @param color Color in a HTML/CSS compatible format, e.g. <strong>#00ff00</strong> or <strong>green</strong> for green.
 * @param opacity Opacity of the fill color. Must be between 0 and 1.
 */
JXG.AbstractRenderer.prototype.setObjectStrokeColor = function(/** JXG.GeometryElement */ obj, /** string */ color, /** number */ opacity) { };

/**
 * Sets an objects fill color.
 * @param obj Reference of the object that wants a new fill color.
 * @param color Color in a HTML/CSS compatible format. If you don't want any fill color at all, choose 'none'.
 * @param opacity Opacity of the fill color. Must be between 0 and 1.
 */
JXG.AbstractRenderer.prototype.setObjectFillColor = function(/** JXG.GeometryElement */ obj, /** string */ color, /** number */ opacity) { };

/**
 * Puts an object into draft mode, i.e. it's visual appearance will be changed. For GEONE<sub>x</sub>T backwards compatibility. 
 * @param obj Reference of the object that shall be in draft mode.
 */
JXG.AbstractRenderer.prototype.setDraft = function (/** JXG.GeometryElement */ obj) {
    if (!obj.visProp['draft']) {
        return;
    }
    var draftColor = obj.board.options.elements.draft.color,
        draftOpacity = obj.board.options.elements.draft.opacity;
        
    if(obj.type == JXG.OBJECTT_TYPE_POLYGON) {
        this.setObjectFillColor(obj, draftColor, draftOpacity);
    }     
    else {
        if(obj.elementClass == JXG.OBJECT_CLASS_POINT) {
            this.setObjectFillColor(obj, draftColor, draftOpacity); 
        }
        else {
            this.setObjectFillColor(obj, 'none', 0); 
        }
        this.setObjectStrokeColor(obj, draftColor, draftOpacity);    
        this.setObjectStrokeWidth(obj, obj.board.options.elements.draft.strokeWidth);
    }      
};

/**
 * Puts an object from draft mode back into normal mode.
 * @param obj Reference of the object that shall no longer be in draft mode.
 */
JXG.AbstractRenderer.prototype.removeDraft = function (/** JXG.GeometryElement */ obj) {
    if(obj.type == JXG.OBJECT_TYPE_POLYGON) {
        this.setObjectFillColor(obj, obj.visProp['fillColor'], obj.visProp['fillColorOpacity']);
    }     
    else {
        if(obj.type == JXG.OBJECT_CLASS_POINT) {
            this.setObjectFillColor(obj, obj.visProp['fillColor'], obj.visProp['fillColorOpacity']);
        }
        this.setObjectStrokeColor(obj, obj.visProp['strokeColor'], obj.visProp['strokeColorOpacity']);        
        this.setObjectStrokeWidth(obj, obj.visProp['strokeWidth']);
    }      
};

/**
 * Highlights an object,
 * i.e. uses the respective highlighting properties of an object.
 * @param obj Reference of the object that will be highlighted.
 */
JXG.AbstractRenderer.prototype.highlight = function(/** JXG.GeometryElement */ obj) {
    var i;
    if(obj.visProp['draft'] == false) {
        if(obj.type == JXG.OBJECT_CLASS_POINT) {
            this.setObjectStrokeColor(obj, obj.visProp['highlightStrokeColor'], obj.visProp['highlightStrokeOpacity']);
            this.setObjectFillColor(obj, obj.visProp['highlightStrokeColor'], obj.visProp['highlightStrokeOpacity']);
        }
        else if(obj.type == JXG.OBJECT_TYPE_POLYGON) {
            this.setObjectFillColor(obj, obj.visProp['highlightFillColor'], obj.visProp['highlightFillOpacity']);
            for(i=0; i<obj.borders.length; i++) {
                this.setObjectStrokeColor(obj.borders[i], obj.borders[i].visProp['highlightStrokeColor'], obj.visProp['highlightStrokeOpacity']);
            }
        }    
        else {
            this.setObjectStrokeColor(obj, obj.visProp['highlightStrokeColor'], obj.visProp['highlightStrokeOpacity']);
            this.setObjectFillColor(obj, obj.visProp['highlightFillColor'], obj.visProp['highlightFillOpacity']);    
        }
    }
};

/**
 * Uses the "normal" colors of an object,
 * i.e. the contrasting function to {@link #highlight}.
 * @param obj Reference of the object that will get its normal colors.
 */
JXG.AbstractRenderer.prototype.noHighlight = function(/** JXG.GeometryElement */ obj) {
    var i;
    if(obj.visProp['draft'] == false) {
        if(obj.type == JXG.OBJECT_CLASS_POINT) {
            this.setObjectStrokeColor(obj, obj.visProp['strokeColor'], obj.visProp['strokeOpacity']);
            this.setObjectFillColor(obj, obj.visProp['strokeColor'], obj.visProp['strokeOpacity']);
        }
        else if(obj.type == JXG.OBJECT_TYPE_POLYGON) {
            this.setObjectFillColor(obj, obj.visProp['fillColor'], obj.visProp['fillOpacity']);
            for(i=0; i<obj.borders.length; i++) {
                this.setObjectStrokeColor(obj.borders[i], obj.borders[i].visProp['strokeColor'], obj.visProp['strokeOpacity']);
            }
        }    
        else {
            this.setObjectStrokeColor(obj, obj.visProp['strokeColor'], obj.visProp['strokeOpacity']);
            this.setObjectFillColor(obj, obj.visProp['fillColor'], obj.visProp['fillOpacity']); 
        }
    }
};

/**
 * Removes an HTML-Element from Canvas. Just a stub.
 * @param node The HTMLElement that shall be removed.
 */
JXG.AbstractRenderer.prototype.remove = function(/** HTMLElement */ node) { };


/* ************************** 
 * general renderer related methods
 * **************************/

/**
 * Stop redraw. This method is called before every update, so a non-vector-graphics based renderer
 * can delete the contents of the drawing panel.
 * @see #unsuspendRedraw
 */
JXG.AbstractRenderer.prototype.suspendRedraw = function() { };

/**
 * Restart redraw. This method is called after updating all the rendering node attributes.
 * @see #suspendRedraw
 */
JXG.AbstractRenderer.prototype.unsuspendRedraw = function() { };

/**
 * The tiny zoom bar shown on the bottom of a board (if {@link JXG.Board#showNavigation} is true).
 * @see #updateText
 */
JXG.AbstractRenderer.prototype.drawZoomBar = function(board) { 
    var doc,
        node,
        node_minus,
        node_100,
        node_plus,
        node_larr,
        node_uarr,
        node_darr,
        node_rarr;

    doc = this.container.ownerDocument;
    node = doc.createElement('div');
    
    //node.setAttribute('id', el.id);
    node.className = 'JXGtext';
    node.style.color = '#aaaaaa';
    node.style.backgroundColor = '#f5f5f5'; 
    node.style.padding = '2px';
    node.style.position = 'absolute';
    node.style.fontSize = '10px';  
    node.style.cursor = 'pointer';
    node.style.zIndex = '100';      
    this.container.appendChild(node);
    node.style.right = '5px'; //(board.canvasWidth-100)+ 'px'; 
    node.style.bottom = '5px'; 
    //node.style.top = (board.canvasHeight-22) + 'px';
    
    node_minus = doc.createElement('span');
    node.appendChild(node_minus);
    node_minus.innerHTML = '&nbsp;&ndash;&nbsp;';
    JXG.addEvent(node_minus, 'click', board.zoomOut, board);

    node_100 = doc.createElement('span');
    node.appendChild(node_100);
    node_100.innerHTML = '&nbsp;o&nbsp;';
    JXG.addEvent(node_100, 'click', board.zoom100, board);
    
    node_plus = doc.createElement('span');
    node.appendChild(node_plus);
    node_plus.innerHTML = '&nbsp;+&nbsp;';
    JXG.addEvent(node_plus, 'click', board.zoomIn, board);
    
    node_larr = doc.createElement('span');
    node.appendChild(node_larr);
    node_larr.innerHTML = '&nbsp;&larr;&nbsp;';
    JXG.addEvent(node_larr, 'click', board.clickLeftArrow, board);
    
    node_uarr = doc.createElement('span');
    node.appendChild(node_uarr);
    node_uarr.innerHTML = '&nbsp;&uarr;&nbsp;';
    JXG.addEvent(node_uarr, 'click', board.clickUpArrow, board);
    
    node_darr = doc.createElement('span');
    node.appendChild(node_darr);
    node_darr.innerHTML = '&nbsp;&darr;&nbsp;';
    JXG.addEvent(node_darr, 'click', board.clickDownArrow, board);
    
    node_rarr = doc.createElement('span');
    node.appendChild(node_rarr);
    node_rarr.innerHTML = '&nbsp;&rarr;&nbsp;';
    JXG.addEvent(node_rarr, 'click', board.clickRightArrow, board);

};

/**
 * Wrapper for getElementById for maybe other renderers which elements are not directly accessible by DOM methods like document.getElementById().
 * @param id Unique identifier for element.
 * @return Reference to an JavaScript object. In case of SVG/VMLRenderer it's a reference to an SVG/VML node.
 */
JXG.AbstractRenderer.prototype.getElementById = function(/** string */ id) /** object */ {
    return document.getElementById(id);
};

/**
 * findSplit() is a subroutine for {@link #RamenDouglasPeuker}.
 * It searches for the point between index i and j which 
 * has the largest distance from the line petween poin_i and point_j.
 **/
JXG.AbstractRenderer.prototype.findSplit = function(pts, i, j) {
    var dist = 0, 
        f = i, 
        d, k, ci, cj, ck,
        x0, y0, x1, y1,
        den, lbda;
        
    if (j-i<2) return [-1.0,0];
    
    ci = pts[i].scrCoords;
    cj = pts[j].scrCoords;
    if (isNaN(ci[1]+ci[2]+cj[1]+cj[2])) return [NaN,j];
    
    for (k=i+1; k<j; k++) {
        ck = pts[k].scrCoords;
        x0 = ck[1]-ci[1];
        y0 = ck[2]-ci[2];
        x1 = cj[1]-ci[1];
        y1 = cj[2]-ci[2];
        den = x1*x1+y1*y1;
        if (den>=JXG.Math.eps) {
            lbda = (x0*x1+y0*y1)/den;
            d = x0*x0+y0*y0 - lbda*(x0*x1+y0*y1);
        } else {
            lbda = 0.0;
            d = x0*x0+y0*y0;
        }
        if (lbda<0.0) {
            d = x0*x0+y0*y0;
        } else if (lbda>1.0) {
            x0 = ck[1]-cj[1];
            y0 = ck[2]-cj[2];
            d = x0*x0+y0*y0;
        }
        if (d>dist) {
            dist = d;
            f = k;
        }
    }
    return [Math.sqrt(dist),f];
};

/**
 * RDB() is a subroutine for {@link #RamenDouglasPeuker}.
 * It runs recursively through the point set and searches the
 * point which has the largest distance from the line between the first point and
 * the last point. If the distance from the line is greater than eps, this point is 
 * included in our new point set otherwise it is discarded.
 * If it is taken, we recursively apply the subroutine to the point set before
 * and after the chosen point.
 */
JXG.AbstractRenderer.prototype.RDP = function(pts,i,j,eps,newPts) {
    var result = this.findSplit(pts,i,j);
    if (result[0]>eps) {
        this.RDP(pts, i,result[1], eps,newPts);
        this.RDP(pts, result[1],j, eps,newPts);
    } else {
        newPts.push(pts[j]);
    }
};

/**
  * Ramen-Douglas-Peuker algorithm.
  * It discards points which are not necessary from the polygonal line defined by the point array
  * pts. The computation is done in screen coordinates.
  * Average runtime is O(nlog(n)), worst case runtime is O(n^2), where n is the number of points.
  */
JXG.AbstractRenderer.prototype.RamenDouglasPeuker = function(pts,eps) {
    var newPts = [], i, k, len;
    //return pts;

    len = pts.length;
    
    // Search for the left most point woithout NaN coordinates
    i = 0;
    while (i<len && isNaN(pts[i].scrCoords[1]+pts[i].scrCoords[2])) {i++;}
    // Search for the right most point woithout NaN coordinates
    k = len - 1;
    while (k>i && isNaN(pts[k].scrCoords[1]+pts[k].scrCoords[2])) {k--;}

    // Exit if nothing is left
    if (i>k || i==len) { return []; }
    
    newPts[0] = pts[i];
    this.RDP(pts,i,k,eps,newPts);
    //this.RDP(pts,0,pts.length-1,eps,newPts);
    //$('debug').innerHTML = newPts.length;
    return newPts;
};

/**
 * Sets the shadow properties to a geometry element.
 * @param {JXG.GeometyElement} element Reference to a geometry object, that should get a shadow
 */
JXG.AbstractRenderer.prototype.setShadow = function(element) {
};


JXG.AbstractRenderer.prototype.updatePathStringPoint = function(el, size, type) {
};

JXG.AbstractRenderer.prototype.eval = function(val) {
    if (typeof val=='function') {
        return val();
    } else {
        return val;
    }
};



/*
    Copyright 2008, 
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

JXG.FileReader = new function() {

this.parseFileContent = function(url, board, format) {
    this.request = false;
    var e;
    try {
        this.request = new XMLHttpRequest();
        if(format.toLowerCase()=='raw')
            this.request.overrideMimeType('text/plain; charset=iso-8859-1');
        else
            this.request.overrideMimeType('text/xml; charset=iso-8859-1');
    } catch (e) {
        try {
            this.request = new ActiveXObject("Msxml2.XMLHTTP");
        } catch (e) {
            try {
                this.request = new ActiveXObject("Microsoft.XMLHTTP");
            } catch (e) {
                this.request = false;
            }
        }
    }
    if (!this.request) {
        alert("AJAX not activated!");
        return;
    }
    this.request.open("GET", url, true);
    if(format.toLowerCase()=='raw') {
        this.cbp = function() {
            var request = this.request;
            if (request.readyState == 4) {
                board(request.responseText);
            }
        }; //).bind(this);        
    } else {
        this.cbp = function() {
            var request = this.request;
            if (request.readyState == 4) {
                this.parseString(request.responseText, board, format, false);
            }
        }; //).bind(this);
    }
    this.cb = JXG.bind(this.cbp,this);
    this.request.onreadystatechange = this.cb;

    try {
        this.request.send(null);
    } catch (e) {
        throw new Error("JSXGraph: problems opening " + url + " !");
    }
}; // end: this.parseFileContent

this.cleanWhitespace = function(el) {
    var cur = el.firstChild;
    while ( cur != null ) {
        if ( cur.nodeType == 3 && ! /\S/.test(cur.nodeValue) ) {
            el.removeChild( cur );
        } else if ( cur.nodeType == 1 ) {
            this.cleanWhitespace( cur );
        }
        cur = cur.nextSibling; 
    }
};

this.stringToXMLTree = function(fileStr) {
    // The string "fileStr" is converted into a XML tree.
    if(typeof DOMParser == "undefined") { 
       // IE workaround, since there is no DOMParser
       DOMParser = function () {};
       DOMParser.prototype.parseFromString = function (str, contentType) {
          if (typeof ActiveXObject != "undefined") {
             var d = new ActiveXObject("MSXML.DomDocument");
             d.loadXML(str);
             return d;
          } 
       };
    }
    var parser=new DOMParser();
    
    var tree = parser.parseFromString(fileStr,"text/xml");
    this.cleanWhitespace(tree);
    return tree;
};

this.parseString = function(fileStr, board, format, isString) {
    // fileStr is a string containing the XML code of the construction
    if (format.toLowerCase()=='geonext') { 
        fileStr = JXG.GeonextReader.prepareString(fileStr);
    }
    if (format.toLowerCase()=='geogebra') {
        // if isString is true, fileStr is a base64 encoded string, otherwise it's the zipped file
    	fileStr = JXG.GeogebraReader.prepareString(fileStr, isString);
    }
    if (format.toLowerCase()=='intergeo') {
    	fileStr = JXG.IntergeoReader.prepareString(fileStr);
    }
    board.xmlString = fileStr;
    var tree = this.stringToXMLTree(fileStr);
    // Now, we can walk through the tree
    this.readElements(tree, board, format);
}; // end this.parse

/**
 * Reading the elements of a geonext or geogebra file
 * @param {} tree expects the content of the parsed geonext file returned by function parseFromString
 * @param {Object} board board object
 */
this.readElements = function(tree, board, format) {
    if (format.toLowerCase()=='geonext') { 
        board.suspendUpdate();
        if(tree.getElementsByTagName('GEONEXT').length != 0) {
            JXG.GeonextReader.readGeonext(tree, board);
        }
        board.unsuspendUpdate();
    }
    else if(tree.getElementsByTagName('geogebra').length != 0) {
        JXG.GeogebraReader.readGeogebra(tree, board);
    }
    else if(format.toLowerCase()=='intergeo') {
         JXG.IntergeoReader.readIntergeo(tree, board);
    }
    board.afterLoad();    
}; // end: this.readElements()

}; // end: FileReader()
/*
    Copyright 2008-2010
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview The Board object is defined in this file. Board controls all properties and methods
 * used to manage a geonext board like adding geometric elements, removing them, managing
 * mouse over, drag & drop of geometric objects etc.
 * @author graphjs
 * @version 0.1
 */

/**
 * Constructs a new Board object.
 * @class This is the Board class. It stores all methods and properties required
 * to manage a geonext board like adding geometric elements, removing them, managing
 * mouse over, drag & drop of geometric objects etc.
 * @constructor
 * @param {String,Object} container The id or reference of the html-element the board is drawn in.
 * @param {JXG.AbstractRenderer} renderer The reference of a geonext renderer.
 * @param {String} id Unique identifier for the board, may be an empty string or null or even undefined.
 * @param {JXG.Coords} origin The coordinates where the origin is placed, in user coordinates.
 * @param {float} zoomX Zoom factor in x-axis direction
 * @param {float} zoomY Zoom factor in y-axis direction
 * @param {int} unitX Units in x-axis direction
 * @param {int} unitY Units in y-axis direction
 * @param {int} canvasWidth  The width of canvas
 * @param {int} canvasHeight The height of canvas
 * @param {bool} showCopyright Display the copyright text
 */
JXG.Board = function(container, renderer, id, origin, zoomX, zoomY, unitX, unitY, canvasWidth, canvasHeight, showCopyright) {
    /**
     * Board is in no special mode, objects are highlighted on mouse over and objects may be
     * clicked to start drag&drop.
     * @type int
     * @private
     * @final
     */
    this.BOARD_MODE_NONE = 0x0000;

    /**
     * Board is in drag mode, objects aren't highlighted on mouse over and the object referenced in
     * drag_obj is updated on mouse movement.
     * @type int
     * @see #drag_obj
     * @private
     * @final
     */
    this.BOARD_MODE_DRAG = 0x0001;

    /**
     * Board is in construction mode, objects are highlighted on mouse over and the behaviour of the board
     * is determined by the construction type stored in the field constructionType.
     * @type int
     * @see #constructionType
     * @private
     * @final
     */
    this.BOARD_MODE_CONSTRUCT = 0x0010;

    /**
     * Board is in move origin mode.
     * @type int
     * @private
     * @final
     */
    this.BOARD_MODE_MOVE_ORIGIN = 0x0002;

    /**
     * Updating is made with low quality, e.g. graphs are evaluated at a lesser amount of points.
     * @type int
     * @see #updateQuality
     * @private
     * @final
     */
    this.BOARD_QUALITY_LOW = 0x1;

    /**
     * Updating is made with high quality, e.g. graphs are evaluated at much more points.
     * @type int
     * @see #updateQuality
     * @private
     * @final
     */
    this.BOARD_QUALITY_HIGH = 0x2;

    /**
     * When the board is in construction mode this construction type says we want to construct a point.
     * @type int
     * @private
     * @final
     */
    this.CONSTRUCTION_TYPE_POINT         = 0x43545054;       // CTPT
    /**
     * When the board is in construction mode this construction type says we want to construct a circle.
     * @type int
     * @private
     * @final
     */
    this.CONSTRUCTION_TYPE_CIRCLE        = 0x4354434C;       // CTCL
    /**
     * When the board is in construction mode this construction type says we want to construct a line.
     * @type int
     * @private
     * @final
     */
    this.CONSTRUCTION_TYPE_LINE          = 0x43544C4E;       // CTLN
    /**
     * When the board is in construction mode this construction type says we want to construct a glider.
     * @type int
     * @private
     * @final
     */
    this.CONSTRUCTION_TYPE_GLIDER        = 0x43544744;       // CTSD
    /**
     * When the board is in construction mode this construction type says we want to construct a midpoint.
     * @type int
     * @private
     * @final
     */
    this.CONSTRUCTION_TYPE_MIDPOINT      = 0x43544D50;       // CTMP
    /**
     * When the board is in construction mode this construction type says we want to construct a perpendicular.
     * @type int
     * @private
     * @final
     */
    this.CONSTRUCTION_TYPE_PERPENDICULAR = 0x43545044;       // CTPD
    /**
     * When the board is in construction mode this construction type says we want to construct a parallel.
     * @type int
     * @private
     * @final
     */
    this.CONSTRUCTION_TYPE_PARALLEL      = 0x4354504C;       // CTPL
    /**
     * When the board is in construction mode this construction type says we want to construct a intersection.
     * @type int
     * @private
     * @final
     */
    this.CONSTRUCTION_TYPE_INTERSECTION  = 0x43544953;       // CTIS

    /**
     * The html-id of the html element containing the board.
     * @type String
     */
    this.container = container;

    /**
     * Pointer to the html element containing the board.
     * @type Object
     */
    this.containerObj = document.getElementById(this.container);
    if (this.containerObj==null) {
        throw new Error("\nJSXGraph: HTML container element '" + (box) + "' not found.");
    }
    //this.containerObj.undoPositioned;  //???

    /**
     * A reference to this boards renderer.
     * @private
     * @type AbstractRenderer
     */
    this.renderer = renderer;

    /**
    * Some standard options
    * @type Options
    */
    //this.options = new JXG.Options();
    this.options = JXG.deepCopy(JXG.Options);
    
    /**
     * Dimension of the board.
     * @private
     * @type int
     */
    this.dimension = 2;

    /**
     * Coordinates of the boards origin
     * @type Coords
     */
    this.origin = {};
    this.origin.usrCoords = [1, 0, 0];
    this.origin.scrCoords = [1, origin[0], origin[1]];

    /**
     * Zoom factor in X direction
     * @type int
     */
    this.zoomX = zoomX;

    /**
     * Zoom factor in Y direction
     * @type int
     */
    this.zoomY = zoomY;

    /**
     * This means the number of pixel which represents
     * one unit in user-coordinates in x direction.
     * @type int
     */
    this.unitX = unitX;

    /**
     * This means the number of pixel which represents
     * one unit in user-coordinates in y direction.
     * @type int
     */
    this.unitY = unitY;

    /**
      * Saving some multiplications
      * @private
      */
    this.stretchX = this.zoomX*this.unitX;
    this.stretchY = this.zoomY*this.unitY;

    /**
     * Canvas Width
     * @type int
     */
    this.canvasWidth = canvasWidth;

    /**
     * Canvas Width
     * @type int
     */
    this.canvasHeight = canvasHeight;

    /**
     * Default font size for labels and texts.
     * @type int
     */
    this.fontSize = this.options.fontSize;

    /**
     * A reference to an object of class Algebra.
     * @see Algebra
     * @private
     * @type Algebra
     */
    this.algebra = new JXG.Algebra(this);

    /* If the given id is not valid, generate an unique id */
    if((id != '') && (id != null) && (typeof document.getElementById(id) != 'undefined'))
        this.id = id;
    else
        this.id = this.generateId();

    /**
     * An array containing all hooked functions.
     * @type Array
     */
    this.hooks = [];

    /**
     * An array containing all other boards that are updated after this board has been updated.
     * @private
     * @type Array
     */
    this.dependentBoards = [];

    /**
     * An associative array containing all geometric objects belonging to the board. Key is the id of the object and value is a reference to the object.
     * @private
     * @type Object
     */
    this.objects = {};

    /**
     * This is used for general purpose animations. Stores all the objects that are currently running an animation.
     * @private
     */
    this.animationObjects = {};

    /**
     * An associative array containing all highlighted geometric objects belonging to the board.
     * @private
     * @type Object
     */
    this.highlightedObjects = {};

    /**
     * Number of objects ever created on this board. This includes every object, even invisible and deleted ones.
     * @private
     * @type int
     */
    this.numObjects = 0;

    /**
     * An associative array to store the objects of the board by name. the name of the object is the key and value is a reference to the object.
     * @type Object
     */
    this.elementsByName = {};

    /**
     * The board mode the board is currently in. Possible values are
     * <ul>
     * <li>Board.BOARD_MODE_NONE</li>
     * <li>Board.BOARD_MODE_DRAG</li>
     * <li>Board.BOARD_MODE_CONSTRUCT</li>
     * </ul>
     * @private
     * @type int
     */
    this.mode = this.BOARD_MODE_NONE;

    /**
     * The update quality of the board. In most cases this is set to Board.BOARD_QUALITY_HIGH when mode is not Board.BOARD_MODE_DRAG
     * and Board.QUALITY_HIGH otherwise. Possible values are
     * <ul>
     * <li>BOARD_QUALITY_LOW</li>
     * <li>BOARD_QUALITY_HIGH</li>
     * </ul>
     * @see #mode
     * @private
     * @type int
     */
    this.updateQuality = this.BOARD_QUALITY_HIGH;

   /**
    * If true updates are skipped
     * @private
    * @type bool
    */
   this.isSuspendedRedraw = false;

   /**
    * The way objects can be dragged. If true, objects can only moved on a predefined grid, if false objects can be moved smoothly almost everywhere.
    * @type bool
    */
   this.snapToGrid = this.options.grid.snapToGrid;

   /**
    * The amount of grid points plus one that fit in one unit of user coordinates in x direction.
    * @type int
    */
   this.gridX = this.options.grid.gridX;

   /**
    * The amount of grid points plus one that fit in one unit of user coordinates in y direction.
    * @type int
    */
   this.gridY = this.options.grid.gridY;

   /**
    * Color of the grid.
    * @type string
    */
   this.gridColor = this.options.grid.gridColor;

   /**
    * Opacity of the grid color, between 0 and 1.
    * @type float
    */
   this.gridOpacity = this.options.grid.gridOpacity;

   /**
    * Determines whether the grid is dashed or not.
    * @type bool
    */
   this.gridDash = this.options.grid.gridDash;

   /**
    * The amount of grid points plus one for snapToGrid that fit in one unit of user coordinates in x direction.
    * @type int
    */
   this.snapSizeX = this.options.grid.snapSizeX;

   /**
    * The amount of grid points plus one for snapToGrid that fit in one unit of user coordinates in y direction.
    * @type int
    */
   this.snapSizeY = this.options.grid.snapSizeY;

   this.calculateSnapSizes();

   /**
    * Visibility of the boards grid.
    * @private
    * @type bool
    */
   this.hasGrid = this.options.grid.hasGrid;

   /**
    * The distance from the mouse to the dragged object in x direction when the user clicked the mouse button.
    * @type int
    * @see drag_dy
    * @see #drag_obj
    * @private
    */
   this.drag_dx = 0;

   /**
    * The distance from the mouse to the dragged object in y direction when the user clicked the mouse button.
    * @type int
    * @see drag_dx
    * @see #drag_obj
    * @private
    */
   this.drag_dy = 0;
   
   /**
     * Absolute position of the mouse pointer in screen pixel from the top left corner
     * of the HTML window.
     */
   this.mousePosAbs = [0,0];

    /**
     * Relative position of the mouse pointer in screen pixel from the top left corner
     * of the JSXGraph canvas (the div element contining the board)-
     */
   this.mousePosRel = [0,0];

   /**
    * A reference to the object that is dragged on the board.
    * @private
    * @type Object
    */
   this.drag_obj = null;

   /**
    * string containing the XML text of the construction.
    * it is set in @see FileReader.parseString.
    * Only useful if a construction from GEONExT, Intergeo, ...
    * is read.
    * @type string
    * @private
    */
   this.xmlString = '';

    /*
    * Display the licence text, @see JSXGraph
    */
    if ( (showCopyright!=null && showCopyright) || (showCopyright==null && this.options.showCopyright) ) {
        this.renderer.displayCopyright(JXG.JSXGraph.licenseText,this.options.fontSize);
    }

   /**
    * Full updates are needed after zoom and axis translates.
    * This saves some time during update
    * @private
    * @type bool
    */
   this.needsFullUpdate = false;

   /**
    * if {reducedUpdate} is set to true, then only the dragged element and few (i.e. 2) following
    * elements are updated during mouse move. On muose up the whole construction is
    * updated. This enables JSXGraph even on very slow devices.
    * @private
    * @type bool
    */
    this.reducedUpdate = false;

   /**
    * If GEONExT constructions are displayed,
    * then this property should be set to true.
    * Then no stdform updates and no dragging
    * of lines, circles and curves is possible.
    * @private
    * @type bool
    */
    this.geonextCompatibilityMode = false;

    if (this.options.text.useASCIIMathML) {
        if (typeof translateASCIIMath != 'undefined') {
            init();
        } else {
            this.options.text.useASCIIMathML = false;
        }
    }
   
   /* Event needs to know which methods to call when mouse is moved or clicked */
   // // Event.observe(this.container, 'mousedown', this.mouseDownListener.bind(this));
   //// Event.observe(this.container, 'mousemove', this.mouseMoveListener.bind(this));
   //Event.observe(document, 'mousedown', this.mouseDownListener.bind(this));
   //Event.observe(this.containerObj, 'mousemove', this.mouseMoveListener.bind(this));

   JXG.addEvent(document,'mousedown', this.mouseDownListener, this);
   JXG.addEvent(this.containerObj, 'mousemove', this.mouseMoveListener, this);

   
   /**
	* iPhone-Events
	*/
	 
   //JXG.addEvent(document,'touchstart', this.touchStartListener, this);
   JXG.addEvent(this.containerObj,'touchstart', this.touchStartListener, this);
   JXG.addEvent(this.containerObj, 'touchmove', this.touchMoveListener, this);
   JXG.addEvent(this.containerObj, 'touchend', this.touchEndListener, this);
};

/**
 * @private
 * Generates unique name for the given object. The result depends on object.type, if the object is a point, just capital characters are used, if it is
 * a line just lower case characters. If object is of type Polygon, lower case prefixed with P_ is used and if it's of type circle, lower case characters
 * prefixed with k_ is used. In any other case, lower case chars prefixed with s_ is used.
 * @param {String,Object} object Reference or id or name of an geometry object that is to be named.
 * @return {String} Unique name for the object.
 */
JXG.Board.prototype.generateName = function(object) {
    if(object.type == JXG.OBJECT_TYPE_TICKS)
        return;

    var possibleNames;
    if(object.elementClass == JXG.OBJECT_CLASS_POINT) {
        // points have capital letters
        possibleNames = ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
                                  'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'];
    }
    else {
        // all other elements get lowercase labels
        possibleNames = ['', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',
                                  'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'];
    }


    // how long the name can be at most
    var maxNameLength = 3;
    var pre = '';
    var nameBase = '';
    var post = '';

    if(object.elementClass == JXG.OBJECT_CLASS_POINT || object.elementClass == JXG.OBJECT_CLASS_LINE) {
    }
    else {
        if(object.type == JXG.OBJECT_TYPE_POLYGON) {
            pre = 'P_{';
            post = '}';
        }
        else if(object.type == JXG.OBJECT_TYPE_CIRCLE) {
            pre = 'k_{';
            post = '}';
        }
        else if(object.type == JXG.OBJECT_TYPE_ANGLE) {
            pre = 'W_{';
            post = '}';
        }
        else {
            pre = 's_{';
            post = '}';
        }
    }
    var indices = [];
    var name = '';
    var tmp = '';

    var i = 0;
    var j = 0;

    for(i=0; i<maxNameLength; i++) {
        indices[i] = 0;
    }

    while (indices[maxNameLength-1] < possibleNames.length) {
        for(indices[0]=1; indices[0]<possibleNames.length; indices[0]++) {
            name = pre;

            for(i=maxNameLength; i>0; i--) {
                name += possibleNames[indices[i-1]];
            }

            if (this.elementsByName[name+post] == null) {
                return name+post;
            }

        }
        indices[0] = possibleNames.length;
        for(i=1; i<maxNameLength; i++) {
            if(indices[i-1] == possibleNames.length) {
                indices[i-1] = 1;
                indices[i]++;
            }
        }
    }

    return '';
};

/**
 * Generates unique id for a board. The result is randomly generated and prefixed with 'gxtBoard'.
 * @return {String} Unique id for a board.
 * @private
 */
JXG.Board.prototype.generateId = function () {
    var r = 1;

    // as long as we don't have an unique id generate a new one
    while(JXG.JSXGraph.boards['gxtBoard' + r] != null) {
        r = Math.round(Math.random()*33);
    }

    return ('gxtBoard' + r);
};

/**
 * Composes the unique id for a board. If the ID is empty ('' or null)
 * a new ID is generated, depending on the object type.
 * Additionally, the id of the label is set.
 * As side effect this.numObjects is updated.
 * @param {Object} object Reference of an geometry object that is to be named.
 * @param {int} object Type of the object.
 * @return {String} Unique id for a board.
 * @private
 */
JXG.Board.prototype.setId = function (obj, type) {
    var num = this.numObjects,
        elId = obj.id;
    this.numObjects++;

    // Falls Id nicht vorgegeben, eine Neue generieren:
    if((elId == '') || (elId == null)) {
        elId = this.id + type + num;
    }
    // Objekt an den Renderer zum Zeichnen uebergeben
    obj.id = elId;
    // Objekt in das assoziative Array einfuegen
    this.objects[elId] = obj;

    if(obj.hasLabel) {
        obj.label.content.id = elId+"Label";
        this.addText(obj.label.content);
    }
    return elId;
};

/**
 * @private
 * Calculates mouse coordinates relative to the boards container.
 * @param {Event} Evt The browsers event object.
 * @type Array
 * @return Array of coordinates relative the boards container top left corner.
 */
JXG.Board.prototype.getRelativeMouseCoordinates = function (Evt) {
    var pCont = this.containerObj,
        cPos = JXG.getOffset(pCont), 
        n; //Element.cumulativeOffset(pCont);

    // add border width
    n = parseInt(JXG.getStyle(pCont,'borderLeftWidth'));
    if (isNaN(n)) n = 0; // IE problem if border-width not set explicitly
    cPos[0] += n;

    n = parseInt(JXG.getStyle(pCont,'borderTopWidth'));
    if (isNaN(n)) n = 0;
    cPos[1] += n;

    // add padding
    n = parseInt(JXG.getStyle(pCont,'paddingLeft'));
    if (isNaN(n)) n = 0;
    cPos[0] += n;

    n = parseInt(JXG.getStyle(pCont,'paddingTop'));
    if (isNaN(n)) n = 0;
    cPos[1] += n;

    return cPos;
};

/**
 * @private
 * Handler for click on left arrow in the navigation bar
 **/
JXG.Board.prototype.clickLeftArrow = function (Event) {
    this.origin.scrCoords[1] += this.canvasWidth*0.1;
    this.moveOrigin();
    return this;
};

/**
 * @private
 * Handler for click on right arrow in the navigation bar
 **/
JXG.Board.prototype.clickRightArrow = function (Event) {
    this.origin.scrCoords[1] -= this.canvasWidth*0.1;
    this.moveOrigin();
    return this;
};

/**
 * @private
 * Handler for click on up arrow in the navigation bar
 **/
JXG.Board.prototype.clickUpArrow = function (Event) {
    this.origin.scrCoords[2] += this.canvasHeight*0.1;
    this.moveOrigin();
    return this;
};

/**
 * @private
 * Handler for click on down arrow in the navigation bar
 **/
JXG.Board.prototype.clickDownArrow = function (Event) {
    this.origin.scrCoords[2] -= this.canvasHeight*0.1;
    this.moveOrigin();
    return this;
};

/**
 * iPhone-Events
 */
 
JXG.Board.prototype.touchStartListener = function (evt) {
	var e = document.createEvent("MouseEvents");   
    this.options.precision.hasPoint = this.options.precision.touch;
	e.initMouseEvent('mousedown', true, false, this.containerObj, 0, evt.targetTouches[0].screenX, evt.targetTouches[0].screenY, evt.targetTouches[0].clientX, evt.targetTouches[0].clientY, false, false, evt.targetTouches.length == 1 ? false: true, false, 0, null);
	this.mouseDownListener(e);
}

JXG.Board.prototype.touchMoveListener = function (evt) {
	evt.preventDefault();	
	var e = document.createEvent("MouseEvents");   
	e.initMouseEvent('mousemove', true, false, this.containerObj, 0, evt.targetTouches[0].screenX, evt.targetTouches[0].screenY, evt.targetTouches[0].clientX, evt.targetTouches[0].clientY, false, false, evt.targetTouches.length == 1 ? false: true, false, 0, null);
	this.mouseMoveListener(e);
}

JXG.Board.prototype.touchEndListener = function (evt) {
	var e = document.createEvent("MouseEvents");   
	e.initMouseEvent('mouseup', true, false, this.containerObj, 0, 0, 0, 0, 0, false, false, false, false, 0, null);
	this.mouseUpListener(e);
    this.options.precision.hasPoint = this.options.precision.mouse;
}

/**
 * This method is called by the browser when the left mouse button is released.
 * @param {Event} Event The browsers event object.
 * @private
 */
JXG.Board.prototype.mouseUpListener = function (evt) {
    // redraw with high precision
    this.updateQuality = this.BOARD_QUALITY_HIGH;

    // release mouseup listener
    JXG.removeEvent(document, 'mouseup', this.mouseUpListener, this);

    this.mode = this.BOARD_MODE_NONE;

    // if origin was moved update everything
    if(this.mode == this.BOARD_MODE_MOVE_ORIGIN) {
        this.moveOrigin();
    } else {
        //this.fullUpdate(); // Full update only needed on moveOrigin? (AW)
        this.update();
    }

    // release dragged object
    this.drag_obj = null;
};

/**
 * This method is called by the browser when the mouse is moved.
 * @param {Event} Evt The browsers event object.
 * @private
 */
JXG.Board.prototype.mouseDownListener = function (Evt) {
    var el, pEl, cPos, absPos, dx, dy;

    cPos = this.getRelativeMouseCoordinates(Evt);
    // position of mouse cursor relative to containers position of container
    absPos = JXG.getPosition(Evt);
    dx = absPos[0]-cPos[0]; //Event.pointerX(Evt) - cPos[0];
    dy = absPos[1]-cPos[1]; //Event.pointerY(Evt) - cPos[1];
    this.mousePosAbs = absPos; // Save the mouse position
    this.mousePosRel = [dx,dy];

    if(Evt.shiftKey) {
        this.drag_dx = dx - this.origin.scrCoords[1];
        this.drag_dy = dy - this.origin.scrCoords[2];
        this.mode = this.BOARD_MODE_MOVE_ORIGIN;
        //Event.observe(this.container, 'mouseup', this.mouseUpListener.bind(this));
        JXG.addEvent(document, 'mouseup', this.mouseUpListener, this);
        return;
    }
    if (this.mode==this.BOARD_MODE_CONSTRUCT) return;

    this.mode = this.BOARD_MODE_DRAG;
    if (this.mode==this.BOARD_MODE_DRAG) {
        for(el in this.objects) {
            pEl = this.objects[el];
            if( (pEl.hasPoint != undefined)
                    && ((pEl.type == JXG.OBJECT_TYPE_POINT) || (pEl.type == JXG.OBJECT_TYPE_GLIDER)
                        /*|| (!this.geonextCompatibilityMode && pEl.type == JXG.OBJECT_TYPE_LINE)  // not yet
                        || (!this.geonextCompatibilityMode && pEl.type == JXG.OBJECT_TYPE_CIRCLE)
                        || (!this.geonextCompatibilityMode && pEl.type == JXG.OBJECT_TYPE_CURVE)*/ )
                    && (pEl.visProp['visible'])
                    && (!pEl.fixed)
                    && (pEl.hasPoint(dx, dy))
                    ) {
                // Points are preferred:
                if ((pEl.type == JXG.OBJECT_TYPE_POINT) || (pEl.type == JXG.OBJECT_TYPE_GLIDER)) {
                    this.drag_obj = this.objects[el];
                    if (this.options.takeFirst) break;
                }
            }
        }
    }

    // if no draggable object can be found, get outta here immediately
    if(this.drag_obj == null) {
        this.mode = this.BOARD_MODE_NONE;
        return;
    }

    /**
      * New mouse position in screen coordinates.
      */
    this.dragObjCoords = new JXG.Coords(JXG.COORDS_BY_SCREEN, [dx,dy], this);
    JXG.addEvent(document, 'mouseup', this.mouseUpListener,this);
};

/**
 * This method is called by the browser when the left mouse button is clicked.
 * @param {Event} Event The browsers event object.
 * @private
 */
JXG.Board.prototype.mouseMoveListener = function (Event) {
    var el, pEl, cPos, absPos, newPos, dx, dy;

    cPos = this.getRelativeMouseCoordinates(Event);
    // position of mouse cursor relative to containers position of container
    absPos = JXG.getPosition(Event);
    dx = absPos[0]-cPos[0]; //Event.pointerX(Evt) - cPos[0];
    dy = absPos[1]-cPos[1]; //Event.pointerY(Evt) - cPos[1];

    this.mousePosAbs = absPos; // Save the mouse position
    this.mousePosRel = [dx,dy];

    this.updateQuality = this.BOARD_QUALITY_LOW;

    this.dehighlightAll(dx,dy);
    if(this.mode != this.BOARD_MODE_DRAG) {
        this.renderer.hide(this.infobox);
    }

    if(this.mode == this.BOARD_MODE_MOVE_ORIGIN) {
        this.origin.scrCoords[1] = dx - this.drag_dx;
        this.origin.scrCoords[2] = dy - this.drag_dy;
        this.moveOrigin();
    }
    else if(this.mode == this.BOARD_MODE_DRAG) {
        newPos = new JXG.Coords(JXG.COORDS_BY_SCREEN, this.getScrCoordsOfMouse(dx,dy), this);
        if (this.drag_obj.type == JXG.OBJECT_TYPE_POINT
            || this.drag_obj.type == JXG.OBJECT_TYPE_LINE
            || this.drag_obj.type == JXG.OBJECT_TYPE_CIRCLE
            || this.drag_obj.type == JXG.OBJECT_TYPE_CURVE) {

/*
            // Do not use setPositionByTransform at the moment!
            // This concept still has to be worked out.
            
            if ((this.geonextCompatibilityMode && this.drag_obj.type==JXG.OBJECT_TYPE_POINT) || this.drag_obj.group.length != 0) {
                // This is for performance reasons with GEONExT files and for groups (transformations do not work yet with groups)
                this.drag_obj.setPositionDirectly(JXG.COORDS_BY_USER,newPos.usrCoords[1],newPos.usrCoords[2]);
            } else {
                this.drag_obj.setPositionByTransform(JXG.COORDS_BY_USER,
                    newPos.usrCoords[1]-this.dragObjCoords.usrCoords[1],
                    newPos.usrCoords[2]-this.dragObjCoords.usrCoords[2]);
                // Save new mouse position in screen coordinates.
                this.dragObjCoords = newPos;
            }
*/            
            this.drag_obj.setPositionDirectly(JXG.COORDS_BY_USER,newPos.usrCoords[1],newPos.usrCoords[2]);
            this.update(this.drag_obj);
        } else if(this.drag_obj.type == JXG.OBJECT_TYPE_GLIDER) {
            var oldCoords = this.drag_obj.coords;
            // First the new position of the glider is set to the new mouse position
            this.drag_obj.setPositionDirectly(JXG.COORDS_BY_USER,newPos.usrCoords[1],newPos.usrCoords[2]);
            // Then, from this position we compute the projection to the object the glider on which the glider lives.
            if(this.drag_obj.slideObject.type == JXG.OBJECT_TYPE_CIRCLE) {
                this.drag_obj.coords = this.algebra.projectPointToCircle(this.drag_obj, this.drag_obj.slideObject);
            } else if (this.drag_obj.slideObject.type == JXG.OBJECT_TYPE_LINE) {
                this.drag_obj.coords = this.algebra.projectPointToLine(this.drag_obj, this.drag_obj.slideObject);
            }
            // Now, we have to adjust the other group elements again.
            if(this.drag_obj.group.length != 0) {
                this.drag_obj.group[this.drag_obj.group.length-1].dX = this.drag_obj.coords.scrCoords[1] - oldCoords.scrCoords[1];
                this.drag_obj.group[this.drag_obj.group.length-1].dY = this.drag_obj.coords.scrCoords[2] - oldCoords.scrCoords[2];
                this.drag_obj.group[this.drag_obj.group.length-1].update(this);
            } else {
                this.update(this.drag_obj);
            }
        }
        this.updateInfobox(this.drag_obj);
    }
    else { // BOARD_MODE_NONE or BOARD_MODE_CONSTRUCT
        // Elements  below the mouse pointer which are not highlighted are highlighted.
        for(el in this.objects) {
            pEl = this.objects[el];
            if( pEl.hasPoint!=undefined && pEl.visProp['visible']==true && pEl.hasPoint(dx, dy)) {
                //this.renderer.highlight(pEl);

                // this is required in any case because otherwise the box won't be shown until the point is dragged
                this.updateInfobox(pEl);
                if(this.highlightedObjects[el] == null) { // highlight only if not highlighted
                    pEl.highlight();
                    this.highlightedObjects[el] = pEl;
                }
            }
        }
    }
    this.updateQuality = this.BOARD_QUALITY_HIGH;
};

/**
 * Updates and displays a little info box to show coordinates of current selected points.
 * @param {JXG.GeometryElement} el A GeometryElement
 * @private
 */
JXG.Board.prototype.updateInfobox = function(el) {
    var x, y, xc, yc;
    if (!el.showInfobox) {
        return this;
    }
    if (el.elementClass == JXG.OBJECT_CLASS_POINT) {
        xc = el.coords.usrCoords[1]*1;
        yc = el.coords.usrCoords[2]*1;
        this.infobox.setCoords(xc+this.infobox.distanceX/(this.stretchX),
                               yc+this.infobox.distanceY/(this.stretchY));
        x = Math.abs(xc);
        if (x>0.1) {
            x = xc.toFixed(2);
        } else if (x>=0.01) {
            x = xc.toFixed(4);
        } else if (x>=0.0001) {
            x = xc.toFixed(6);
        } else {
            x = xc;
        }
        y = Math.abs(yc);
        if (y>0.1) {
            y = yc.toFixed(2);
        } else if (y>=0.01) {
            y = yc.toFixed(4);
        } else if (y>=0.0001) {
            y = yc.toFixed(6);
        } else {
            y = yc;
        }

        this.highlightInfobox(x,y,el);
        this.renderer.show(this.infobox);
        this.renderer.updateText(this.infobox);
    }
    return this;
};

JXG.Board.prototype.highlightInfobox = function(x,y,el) {
    this.infobox.setText('<span style="color:#bbbbbb;">(' + x + ', ' + y + ')</span>');
    return this;
};

/**
 * Remove highlighting of all elements.
 * @private
 */
JXG.Board.prototype.dehighlightAll = function(x,y) {
    var el, pEl;

    for(el in this.highlightedObjects) {
        //this.renderer.noHighlight(this.highlightedObjects[el]);
        pEl = this.highlightedObjects[el];
        if((pEl.hasPoint == undefined) ||
           (!pEl.hasPoint(x, y)) ||
           (pEl.visProp['visible'] == false)) { // dehighlight only if necessary
                pEl.noHighlight();
                delete(this.highlightedObjects[el]);
        }
    }
    return this;
};

/**
 * In case of snapToGrid activated this method caclulates the screen coords of mouse "snapped to grid".
 * @param {int} x X coordinate in screen coordinates
 * @param {int} y Y coordinate in screen coordinates
 */
JXG.Board.prototype.getScrCoordsOfMouse = function (x,y) {
    if(this.snapToGrid) {
        var newCoords = new JXG.Coords(JXG.COORDS_BY_SCREEN, [x,y], this);
        newCoords.setCoordinates(JXG.COORDS_BY_USER,
            [Math.round((newCoords.usrCoords[1])*this.snapSizeX)/this.snapSizeX,
             Math.round((newCoords.usrCoords[2])*this.snapSizeY)/this.snapSizeY]);
        return [newCoords.scrCoords[1], newCoords.scrCoords[2]];
    } else {
        return [x,y];
    }
};

/**
 * In case of snapToGrid activated this method caclulates the user coords of mouse "snapped to grid".
 * @param {Event} Evt Event object containing the mouse coordinates.
 */
JXG.Board.prototype.getUsrCoordsOfMouse = function (Evt) {
    var cPos = this.getRelativeMouseCoordinates(Evt);
    //var x = Event.pointerX(Evt) - cPos[0];
    //var y = Event.pointerY(Evt) - cPos[1];
    var absPos = JXG.getPosition(Evt);
    var x = absPos[0]-cPos[0]; //Event.pointerX(Evt) - cPos[0];
    var y = absPos[1]-cPos[1]; //Event.pointerY(Evt) - cPos[1];

    var newCoords = new JXG.Coords(JXG.COORDS_BY_SCREEN, [x,y], this);
    if(this.snapToGrid) {
        newCoords.setCoordinates(JXG.COORDS_BY_USER,
            [Math.round((newCoords.usrCoords[1])*this.snapSizeX)/this.snapSizeX,
             Math.round((newCoords.usrCoords[2])*this.snapSizeY)/this.snapSizeY]);
    }
    return [newCoords.usrCoords[1], newCoords.usrCoords[2]];
};

/**
 * Collects all elements under current mouse position plus current user coordinates of mouse cursor.
 * @param {Event} Evt Event object containing the mouse coordinates.
 * @type Array
 * @return Array of elements at the current mouse position plus current user coordinates of mouse.
 * @private
 */
JXG.Board.prototype.getAllUnderMouse = function (Evt) {
    var elList = this.getAllObjectsUnderMouse(Evt);
    elList.push(this.getUsrCoordsOfMouse(Evt));
    return elList;
    //return {"elList":elList, "coords":this.getUsrCoordsOfMouse(Evt)};
};

/**
 * Collects all elements under current mouse position.
 * @param {Event} Evt Event object containing the mouse coordinates.
 * @type Array
 * @return Array of elements at the current mouse position.
 * @private
 */
JXG.Board.prototype.getAllObjectsUnderMouse = function (Evt) {
    var cPos = this.getRelativeMouseCoordinates(Evt);

    // mouse position relative to container
    //var dx = Event.pointerX(Evt) - cPos[0];
    //var dy = Event.pointerY(Evt) - cPos[1];
    var absPos = JXG.getPosition(Evt);
    var dx = absPos[0]-cPos[0]; //Event.pointerX(Evt) - cPos[0];
    var dy = absPos[1]-cPos[1]; //Event.pointerY(Evt) - cPos[1];
    var elList = [];
    for (var el in this.objects) {
        if (this.objects[el].visProp['visible'] && this.objects[el].hasPoint(dx, dy)) {
            elList.push(this.objects[el]);
        }
    }
    return elList;
};

/**
 * Sets the board mode.
 * @param {int} mode The board mode the board should be set to. Possible values are
 * <li><ul>BOARD_MODE_NONE</ul><ul>BOARD_MODE_DRAG</ul><ul>BOARD_MODE_CONSTRUCT</ul><ul>BOARD_MODE_MOVE_ORIGIN</ul></li>
 * @private
 */
JXG.Board.prototype.setBoardMode = function (mode) {
    this.mode = mode;
    return this;
};

/**
 * Moves the origin and initializes an update of all elements.
 * @private
 */
JXG.Board.prototype.moveOrigin = function () {
    for(var Element in this.objects) {
        if( (this.objects[Element].elementClass == JXG.OBJECT_CLASS_POINT) ||
            (this.objects[Element].type == JXG.OBJECT_TYPE_CURVE) ||
            (this.objects[Element].type == JXG.OBJECT_TYPE_AXIS) ||
            (this.objects[Element].type == JXG.OBJECT_TYPE_TEXT) ) {
            if((this.objects[Element].type != JXG.OBJECT_TYPE_CURVE) && (this.objects[Element].type != JXG.OBJECT_TYPE_AXIS))
                this.objects[Element].coords.usr2screen();
        }
    }

    this.clearTraces();

    this.fullUpdate();
    if(this.hasGrid) {
        this.renderer.removeGrid(this);
        this.renderer.drawGrid(this);
    }
    return this;
};

/**
 * After construction of the object the visibility is set
 * and the label is constructed if necessary.
 * @param {Object} obj The object to add.
 * @private
 */
JXG.Board.prototype.finalizeAdding = function (obj) {
    if (obj.hasLabel) {
        this.renderer.drawText(obj.label.content);
    }
    if(!obj.visProp['visible']) {
        this.renderer.hide(obj);
    }

    if(obj.hasLabel && !obj.label.content.visProp['visible']) {
        this.renderer.hide(obj.label.content);
    }
};

/**
 * Registers a point at the board and adds it to the renderer.
 * @param {JXG.Point} obj The point to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addPoint = function (obj) {
    //this.elementsByName[obj.name] = obj;
    var id = this.setId(obj,'P');
    this.renderer.drawPoint(obj);
    this.finalizeAdding(obj);
    return id;
};

/**
 * Registers a line at the board and adds it to the renderer.
 * @param {JXG.Line} obj The line to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addLine = function (obj) {
    var id = this.setId(obj,'L');
    this.renderer.drawLine(obj);
    this.finalizeAdding(obj);
    return id;
};

/**
 * Registers a circle at the board and adds it to the renderer.
 * @param {JXG.Circle} obj The circle to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addCircle = function(obj) {
    var id = this.setId(obj,'C');
    this.renderer.drawCircle(obj);
    this.finalizeAdding(obj);
    return id;
};

/**
 * Registers a polygon at the board and adds it to the renderer.
 * @param {JXG.Polygon} obj The polygon to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addPolygon = function(obj) {
    var id = this.setId(obj,'Py');
    this.renderer.drawPolygon(obj);
    this.finalizeAdding(obj);
    return id;
};

/**
 * Registers a arc at the board and adds it to the renderer.
 * @param {JXG.Arc} obj The arc to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addArc = function(obj) {
    var id = this.setId(obj,'Ac');
    this.renderer.drawArc(obj);
    this.finalizeAdding(obj);
    return id;
};

/**
 * Registers a sector at the board and adds it to the renderer.
 * @param {JXG.Sector} obj The sector to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addSector = function(obj) {
    return this.setId(obj,'Sc');
};

/**
 * Registers an angle at the board and adds it to the renderer.
 * @param {JXG.Angle} obj The angle to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addAngle = function (obj) {
    var id = this.setId(obj,'Ag');
    this.renderer.drawAngle(obj);
    this.finalizeAdding(obj);
    return id;
};

/**
 * Registers a curve at the board and adds it to the renderer.
 * @param {JXG.Curve} obj The curve to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addCurve = function (obj) {
    var id = this.setId(obj,'G');
    this.renderer.drawCurve(obj);
    this.finalizeAdding(obj);
    return id;
};

/**
 * Registers a chart at the board and adds it to the renderer.
 * @param {JXG.Chart} obj The chart to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addChart = function (obj) {
    return this.setId(obj,'Chart');
};

/**
 * Registers a arrow at the board and adds it to the renderer.
 * @param {JXG.Arrow} obj The arrow to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addArrow = function(obj) {
    var num = this.numObjects, elId;
    this.numObjects++;

    // Falls Id nicht vorgegeben, eine Neue generieren:
    elId = obj.id;
    if((elId == '') || (elId == null)) {
        elId = this.id + 'A' + num;
    }

    // Objekt in das assoziative Array einfuegen
    this.objects[elId] = obj;

    // Objekt an den Renderer zum Zeichnen uebergeben
    obj.id = elId;
    this.renderer.drawArrow(obj);

    return elId;
};

/**
 * Adds a line to the board and renderer which is orthogonal to the given line and contains point.
 * @param {JXG.Line} l A line.
 * @param {JXG.Point} p A Point.
 * @param {String} id Unique identifier for this object.  If null or an empty string is given,
 * an unique id will be generated by Board
 * @param {String} name Not necessarily unique name, displayed on the board.  If null or an
 * empty string is given, an unique name will be generated.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addNormal = function(l, p, id, name) {
    var point = JXG.getReference(this, p);
    var line = JXG.getReference(this, l);

    var number = this.numObjects;
    number++;
    if((id == '') || (id == null)) {
        id = this.id + 'L' + number;
    }

    // versteckter Hilfs-Punkt
    var erg = this.algebra.perpendicular(line, point);
    var p2coords = erg[0].usrCoords.slice(1);
    var point2 = new JXG.Point(this, p2coords, id+"P2", '', false);
    point2.fixed = true;
    point.addChild(point2); // notwendig, um auch den Punkt upzudaten
    line.addChild(point2); // notwendig, um auch den Punkt upzudaten

    var perpendicular;
    if(erg[1]) {
        perpendicular = new JXG.Line(this, point2.id, point.id, id, name);
    }
    else {
        perpendicular = new JXG.Line(this, point.id, point2.id, id, name);
    }
    perpendicular.changed = erg[1];
    //point.addChild(perpendicular);
    //line.addChild(perpendicular);

    perpendicular.update = function() {
        if (this.needsUpdate) {
            var erg = this.board.algebra.perpendicular(line, point);
            point2.coords = erg[0];
            if(this.changed != erg[1]) {
                var tmp = this.point1;
                this.point1 = this.point2;
                this.point2 = tmp;
            }
            this.updateStdform(); // For the new intersection functions
            if(this.traced) {
                this.cloneToBackground(true);
            }
        }
    };
    return perpendicular;
};

/**
 * Registers an intersection at the board and adds it to the renderer.
 * @param {JXG.Intersection} obj The intersection to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addIntersection = function (obj) {
    var number = this.numObjects;
    this.numObjects++;
    var elementId = obj.id;

    // Falls Id nicht vergeben, eine neue generieren:
    if((elementId == '') || (elementId == null)) {
        elementId = this.id + 'I' + number;
    }

    // Objekt in das assoziative Array einfuegen
    this.objects[elementId] = obj;

    obj.id = elementId;

    obj.intersect1.addChild(obj);
    obj.intersect2.addChild(obj);

    return elementId;
};

/**
 * Registers a text at the board and adds it to the renderer.
 * @param {JXG.Text} obj The text to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addText = function (obj) {
    var number = this.numObjects;
    this.numObjects++;

    // Falls Id nicht vergeben, eine Neue generieren:
    var elementId = obj.id;
    if((elementId == '') || (elementId == null)) {
        elementId = this.id + 'T' + number;
    }

    // Objekt in das assoziative Array einfuegen
    this.objects[elementId] = obj;

    // Objekt an den Renderer zum Zeichnen uebergeben
    obj.id = elementId;
    if(!obj.isLabel) {
        this.renderer.drawText(obj);
        if(!obj.visProp['visible']) {
            this.renderer.hide(obj);
        }
    }

    return elementId;
};

/**
  * Add conditional updates to the elements.
  * @param {string} str String containing coniditional update in geonext syntax
  */
JXG.Board.prototype.addConditions = function (str) {
    var res = null;
    var plaintext = 'var el,x,y,c;\n';
    var i = str.indexOf('<data>');
    var j = str.indexOf('</data>');
    if (i<0) {
        return;
    }
    while (i>=0) {
        var term = str.slice(i+6,j); // throw away <data>
        var m = term.indexOf('=');
        var left = term.slice(0,m);
        var right = term.slice(m+1);
        m = left.indexOf('.'); // Dies erzeugt Probleme bei Variablennamen der Form " Steuern akt."
        var name = left.slice(0,m);    //.replace(/\s+$/,''); // do NOT cut out name (with whitespace)
        var el = this.elementsByName[JXG.unescapeHTML(name)];

        var property = left.slice(m+1).replace(/\s+/g,'').toLowerCase(); // remove whitespace in property
        right = this.algebra.geonext2JS(right);
        right = right.replace(/this\.board\./g,'this.');

        // Debug
        if (typeof this.elementsByName[name]=='undefined'){
            alert("debug conditions: |"+name+"| undefined");
        }
        plaintext += "el = this.objects[\"" + el.id + "\"];\n";
        //plaintext += "if (el==undefined) { $('debug').value = \"" + name + "\"; } else {\n";
        switch (property) {
            case 'x':
                plaintext += 'y=el.coords.usrCoords[2];\n';  // y stays
                //plaintext += 'el.coords=new JXG.Coords(JXG.COORDS_BY_USER,['+(right) +',y],this);\n';
                plaintext += 'el.setPositionDirectly(JXG.COORDS_BY_USER,'+(right) +',y);\n';
                plaintext += 'el.update();\n';
                break;
            case 'y':
                plaintext += 'x=el.coords.usrCoords[1];\n';  // x stays
                plaintext += 'el.coords=new JXG.Coords(JXG.COORDS_BY_USER,[x,'+(right)+'],this);\n';
                //plaintext += 'el.update();\n';
                break;
            case 'visible':
                plaintext += 'c='+(right)+';\n';
                plaintext += 'if (c) {el.showElement();} else {el.hideElement();}\n';
                break;
            case 'position':
                plaintext += 'el.position = ' + (right) +';\n';
                plaintext += 'el.update();\n';
                //plaintext += 'this.updateElements();\n';
                break;
            case 'stroke':
                plaintext += 'el.strokeColor = ' + (right) +';\n';
                break;
            case 'strokewidth':
                plaintext += 'el.strokeWidth = ' + (right) +';\n';   // wird auch bei Punkten verwendet, was nicht realisiert ist.
                break;
            case 'label':
                //plaintext += 'var color = ' + (right) +';\n';
                //plaintext += 'el.setProperty("labelColor:color");\n';
                break;
            default:
                alert("property '" + property + "' in conditions not implemented:" + right);
                break;
        }
        //plaintext += "}\n";
        str = str.slice(j+7); // cut off "</data>"
        i = str.indexOf('<data>');
        j = str.indexOf('</data>');
    }
    plaintext += 'this.prepareUpdate();\n';
    plaintext += 'this.updateElements();\n';
    plaintext += 'return true;\n';
    //alert(plaintext);
    this.updateConditions = new Function(plaintext);
    this.updateConditions();
};

/**
 * Computes the commands in the conditions-section of the gxt file.
 * It is evaluated after an update, before the unsuspendRedraw.
 * The function is generated in @see #addConditions
 * @private
 */
JXG.Board.prototype.updateConditions = function() { return false; };

/**
 * Registers an image at the board and adds it to the renderer.
 * @param {JXG.Image} obj The image to add.
 * @type String
 * @return Element id of the object.
 * @private
 */
JXG.Board.prototype.addImage = function (obj) {
    var number = this.numObjects;
    this.numObjects++;
    var elementId = obj.id;

    // Falls Id nicht vergeben, eine neue generieren:
    if((elementId == '') || (elementId == null)) {
        elementId = this.id + 'Im' + number;
    }

    // Objekt in die assoziativen Arrays einfuegen
    this.objects[elementId] = obj;
    this.elementsByName[obj.name] = obj;

    // Objekt an den Renderer zum Zeichnen uebergeben
    obj.id = elementId;

    this.renderer.drawImage(obj);
    if(!obj.visProp['visible']) {
       this.renderer.hide(obj);
    }

    return elementId;
};

/**
 * Calculates adequate snap sizes.
 * @private
 */
JXG.Board.prototype.calculateSnapSizes = function() {
    var p1 = new JXG.Coords(JXG.COORDS_BY_USER,[0,0],this),
        p2 = new JXG.Coords(JXG.COORDS_BY_USER,[1/this.gridX,1/this.gridY],this),
        x = p1.scrCoords[1]-p2.scrCoords[1],
        y = p1.scrCoords[2]-p2.scrCoords[2];

    this.snapSizeX = this.gridX;
    while(Math.abs(x) > 25) {
        this.snapSizeX *= 2;
        x /= 2;
    }

    this.snapSizeY = this.gridY;
    while(Math.abs(y) > 25) {
        this.snapSizeY *= 2;
        y /= 2;
    }
    return this;
};

/**
 * Apply update on all objects with the
 * new zoom-factors.
 * @private
 */
JXG.Board.prototype.applyZoom = function() {
    var el;

    for(el in this.objects) {
        if( (this.objects[el].elementClass == JXG.OBJECT_CLASS_POINT) ||
            (this.objects[el].type == JXG.OBJECT_TYPE_CURVE) ||
            (this.objects[el].type == JXG.OBJECT_TYPE_AXIS) ||
            (this.objects[el].type == JXG.OBJECT_TYPE_TEXT) ) {
            if((this.objects[el].type != JXG.OBJECT_TYPE_CURVE) && (this.objects[el].type != JXG.OBJECT_TYPE_AXIS))
                this.objects[el].coords.usr2screen();
        }
    }
    this.calculateSnapSizes();

    this.clearTraces();

    this.fullUpdate();
    if(this.hasGrid) {
        this.renderer.removeGrid(this);
        this.renderer.drawGrid(this);
    }
    return this;
};

/**
 * Zooms into the board.
 */
JXG.Board.prototype.zoomIn = function() {
    var oX, oY;
    this.zoomX *= this.options.zoom.factor;
    this.zoomY *= this.options.zoom.factor;
    oX = this.origin.scrCoords[1]*this.options.zoom.factor;
    oY = this.origin.scrCoords[2]*this.options.zoom.factor;
    this.origin = new JXG.Coords(JXG.COORDS_BY_SCREEN, [oX, oY], this);
    this.stretchX = this.zoomX*this.unitX;
    this.stretchY = this.zoomY*this.unitY;
    this.applyZoom();
    return this;
};

/**
 * Zooms out of the board.
 */
JXG.Board.prototype.zoomOut = function() {
    var oX, oY;
    this.zoomX /= this.options.zoom.factor;
    this.zoomY /= this.options.zoom.factor;
    oX = this.origin.scrCoords[1]/this.options.zoom.factor;
    oY = this.origin.scrCoords[2]/this.options.zoom.factor;
    this.origin = new JXG.Coords(JXG.COORDS_BY_SCREEN, [oX, oY], this);
    
    this.stretchX = this.zoomX*this.unitX;
    this.stretchY = this.zoomY*this.unitY;
    this.applyZoom();
    return this;
};

/**
 * Resets zoom factor zu 1.
 */
JXG.Board.prototype.zoom100 = function() {
    var oX, oY, zX, zY;
    
    zX = this.zoomX;
    zY = this.zoomY;
    this.zoomX = 1.0;
    this.zoomY = 1.0;

    oX = this.origin.scrCoords[1]/zX;
    oY = this.origin.scrCoords[2]/zY;
    this.origin = new JXG.Coords(JXG.COORDS_BY_SCREEN, [oX, oY], this);

    this.stretchX = this.zoomX*this.unitX;
    this.stretchY = this.zoomY*this.unitY;
    this.applyZoom();
    return this;
};

/**
 * Zooms the board so every visible point is shown. Keeps aspect ratio.
 */
JXG.Board.prototype.zoomAllPoints = function() {
    var ratio, minX, maxX, minY, maxY, el,
        border, borderX, borderY, distX, distY, newZoom, newZoomX, newZoomY,
        newOriginX, newOriginY;

    ratio = this.zoomX / this.zoomY;
    minX = 0; // (0,0) soll auch sichtbar bleiben
    maxX = 0;
    minY = 0;
    maxY = 0;
    for(el in this.objects) {
        if( (this.objects[el].elementClass == JXG.OBJECT_CLASS_POINT) &&
            this.objects[el].visProp['visible']) {
            if(this.objects[el].coords.usrCoords[1] < minX) {
                minX = this.objects[el].coords.usrCoords[1];
            } else if(this.objects[el].coords.usrCoords[1] > maxX) {
                maxX = this.objects[el].coords.usrCoords[1];
            }
            if(this.objects[el].coords.usrCoords[2] > maxY) {
                maxY = this.objects[el].coords.usrCoords[2];
            } else if(this.objects[el].coords.usrCoords[2] < minY) {
                minY = this.objects[el].coords.usrCoords[2];
            }
        }
    }
    border = 50;
    borderX = border/(this.unitX*this.zoomX);
    borderY = border/(this.unitY*this.zoomY);

    distX = maxX - minX + 2*borderX;
    distY = maxY - minY + 2*borderY;

    newZoom = Math.min(this.canvasWidth/(this.unitX*distX), this.canvasHeight/(this.unitY*distY));
    newZoomY = newZoom;
    newZoomX = newZoom*ratio;

    newOriginX = -(minX-borderX)*this.unitX*newZoomX;
    newOriginY = (maxY+borderY)*this.unitY*newZoomY;
    this.origin = new JXG.Coords(JXG.COORDS_BY_SCREEN, [newOriginX, newOriginY], this);
    this.zoomX = newZoomX;
    this.zoomY = newZoomY;
    this.stretchX = this.zoomX*this.unitX;
    this.stretchY = this.zoomY*this.unitY;

    this.applyZoom();
    return this;
};

/**
 * Removes object from board and renderer.
 * @param {GeometryElement} object The object to remove.
 */
JXG.Board.prototype.removeObject = function(object) {
    var el, i;

    if(JXG.isArray(object)) {
        for(i=0; i<object.length; i++)
            this.removeObject(object[i]);
    }

    object = JXG.getReference(this, object);

    /* Wenn weder die ID noch der Name des Objekts bekannt ist, einfach wieder zurueckgehen */
    if(object == undefined) {
        return this;
    }

    try{
        /* Alle Kinder entfernen */
        for(el in object.childElements) {
            object.childElements[el].board.removeObject(object.childElements[el]);
        }

        for(el in this.objects) {
            if(typeof this.objects[el].childElements != 'undefined')
                delete(this.objects[el].childElements[object.id]);
        }

        /* Das Objekt selbst aus board.objects und board.elementsByName loeschen */
        delete(this.objects[object.id]);
        delete(this.elementsByName[object.name]);

        /* Alles weitere erledigt das Objekt selbst fuer uns. Ist sinnvoller, weil man sonst wieder unterscheiden muesste, was das fuer ein Objekt ist. */
        if(object.remove != undefined) object.remove();
    } catch(e) {
//        alert(object.id + ': Could not be removed, JS says:\n\n' + e);
    }
    return this;
};

/**
 * Initialise some objects which are contained in every GEONExT construction by default,
 * but are not contained in the gxt files.
 * @private
 */
JXG.Board.prototype.initGeonextBoard = function() {
    var p1, p2, p3, l1, l2;

    p1 = new JXG.Point(this, [0,0],this.id + 'gOOe0','Ursprung',false);
    p1.fixed = true;
    p2 = new JXG.Point(this, [1,0],this.id + 'gXOe0','Punkt_1_0',false);
    p2.fixed = true;
    p3 = new JXG.Point(this, [0,1],this.id + 'gYOe0','Punkt_0_1',false);
    p3.fixed = true;
    l1 = new JXG.Line(this, this.id + 'gOOe0', this.id + 'gXOe0', this.id + 'gXLe0','X-Achse');
    l1.hideElement();
    l2 = new JXG.Line(this, this.id + 'gOOe0', this.id + 'gYOe0', this.id + 'gYLe0','Y-Achse');
    l2.hideElement();
    return this;
};

/**
 * Initialise the info box object which is used to display
 * the coordinates of points under the mouse pointer,
 * @private
 */
JXG.Board.prototype.initInfobox= function() {
    //this.infobox = new JXG.Label(this, '0,0', new JXG.Coords(JXG.COORDS_BY_USER, [0, 0], this), this.id + '__infobox');
    this.infobox = new JXG.Text(this, '0,0', '', [0,0], this.id + '__infobox',null, null, false, 'html');
    this.infobox.distanceX = -20;
    this.infobox.distanceY = 25;
    //this.renderer.drawText(this.infobox);
    this.renderer.hide(this.infobox);
    return this;
};

/**
 * Change the height and width of the board's container.
 * @param {int} canvasWidth New width of the container.
 * @param {int} canvasHeight New height of the container.
 */
JXG.Board.prototype.resizeContainer = function(canvasWidth, canvasHeight) {
    this.canvasWidth = 1*canvasWidth;
    this.canvasHeight = 1*canvasHeight;
    this.containerObj.style.width = (this.canvasWidth) + 'px';
    this.containerObj.style.height = (this.canvasHeight) + 'px';
    return this;
};

/**
 * Lists the dependencies graph in a new HTML-window.
 */
JXG.Board.prototype.showDependencies = function() {
    var el, t, c, f, i;

    t = '<p>\n';
    for (el in this.objects) {
        i = 0;
        for (c in this.objects[el].childElements) {
            i++;
        }
        if (i>=0) {
            t += '<b>' + this.objects[el].id + ':</b> ';
        }
        for (c in this.objects[el].childElements) {
            t += this.objects[el].childElements[c].id+'('+this.objects[el].childElements[c].name+')'+', ';
        }
        t += '<p>\n';
    }
    t += '</p>\n';
    f = window.open();
    f.document.open();
    f.document.write(t);
    f.document.close();
    return this;
};

/**
 * Lists the XML code of the construction in a new HTML-window.
 */
JXG.Board.prototype.showXML = function() {
    var f = window.open("");
    f.document.open();
    f.document.write("<pre>"+JXG.escapeHTML(this.xmlString)+"</pre>");
    f.document.close();
    return this;
};

/**
 * Sets for all objects the needsUpdate flag to "true".
 * @param {Object,String} drag Element that caused the update.
 * @private
 */
JXG.Board.prototype.prepareUpdate = function(drag) {
    var el;
    for(el in this.objects) {
       this.objects[el].needsUpdate = true;
    }
    return this;
};

/**
  * Runs through all elements and calls their update() method.
  * @param {Object,String} drag Element that caused the update.
  * @private
  */
JXG.Board.prototype.updateElements = function(drag) {
    var el, pEl,
        isBeforeDrag = true; // If possible, we start the update at the dragged object.

    drag = JXG.getReference(this, drag);
    if (drag==null) {
        isBeforeDrag = false;
    }

    for(el in this.objects) {
        pEl = this.objects[el];
        if (drag!=null && pEl.id != drag.id) {
            isBeforeDrag = false;
        }
        if (!(isBeforeDrag || this.needsFullUpdate || pEl.needsRegularUpdate)) { continue; }
        if (drag==null || pEl.id!=drag.id) {
            //if (this.needsFullUpdate) { pEl.update(true); }
            pEl.update(true);
        } else {
            pEl.update(false);
        }
    }
    return this;
};

/**
  * Runs through all elements and calls their update() method.
  * @param {Object,String} drag Element that caused the update.
  * @private
  */
JXG.Board.prototype.updateRenderer = function(drag) {
    var el, pEl;
    drag = JXG.getReference(this, drag);
    for(el in this.objects) {
        pEl = this.objects[el];
        if (!this.needsFullUpdate && !pEl.needsRegularUpdate) { continue; }
        if (drag == null || pEl.id != drag.id) {
            //if (this.needsFullUpdate) { pEl.updateRenderer(); }
            pEl.updateRenderer();
        } else {
            pEl.updateRenderer();
        }
    }
    return this;
};

/**
  * Adds a hook to this board.
  * @param {function} hook A function to be called by the board after an update occured.
  * @type int
  * @return Id of the hook, required to remove the hook from the board.
  */
JXG.Board.prototype.addHook = function(hook) {
    this.hooks.push(hook);

    hook(this);

    return (this.hooks.length-1);
};

/**
  * Deletes a hook from this board.
  * @param {int} id Id for the hook, required to delete the hook.
  */
JXG.Board.prototype.removeHook = function(id) {
    this.hooks[id] = null;
    return this;
};

/**
  * Runs through all hooked functions and calls them.
  * @private
  */
JXG.Board.prototype.updateHooks = function() {
    var i;
    for(i=0; i<this.hooks.length; i++) {
        if(this.hooks[i] != null)
            this.hooks[i](this);
    }
    return this;
};

/**
  * Adds a dependent board to this board.
  * @param {object}  A reference to board which will be updated after an update of this board occured.
  */
JXG.Board.prototype.addChild = function(board) {
    this.dependentBoards.push(board);
    this.update();
    return this;
};

/**
  * Deletes a board from the list of dependent boards.
  * @param {object} board Reference to the board which will be removed.
  */
JXG.Board.prototype.removeChild = function(board) {
    var i;
    for (i=this.dependentBoards.length-1; i>=0; i--) {
        if (this.dependentBoards[i] == board) {
            this.dependentBoards.splice(i,1);
        }
    }
    return this;
};

/**
  * Runs through most elements and calls their
  * update() method and update the conditions.
  * @param {Object,String} drag Element that caused the update.
  */
JXG.Board.prototype.update = function(drag) {
    var i, len, boardId;

    if (this.isSuspendedUpdate) { return this; }
    this.prepareUpdate(drag).updateElements(drag).updateConditions();
    this.renderer.suspendRedraw();
    this.updateRenderer(drag);
    this.renderer.unsuspendRedraw();
    this.updateHooks();

    // To resolve dependencies between boards
    //for(var board in JXG.JSXGraph.boards) {
    len = this.dependentBoards.length;
    for (i=0; i<len; i++) {
        boardId = this.dependentBoards[i].id;
        if(JXG.JSXGraph.boards[boardId] != this) {
            JXG.JSXGraph.boards[boardId].updateQuality = this.updateQuality;
            JXG.JSXGraph.boards[boardId].prepareUpdate(drag).updateElements(drag).updateConditions();
            JXG.JSXGraph.boards[boardId].renderer.suspendRedraw();
            JXG.JSXGraph.boards[boardId].updateRenderer(drag);
            JXG.JSXGraph.boards[boardId].renderer.unsuspendRedraw();
            JXG.JSXGraph.boards[boardId].updateHooks();
        }

    }
    return this;
};

/**
  * Runs through all elements and calls their
  * update() method and update the conditions.
  * This is necessary after zooming and changing the bounding box.
  */
JXG.Board.prototype.fullUpdate = function() {
    this.needsFullUpdate = true;
    this.update();
    this.needsFullUpdate = false;
    return this;
};

/**
 * Creates a new geometric element of type elementType.
 * @param {string} elementType Type of the element to be constructed given as a string e.g. 'point' or 'circle'.
 * @param {Array} parents Array of parent elements needed to construct the element e.g. coordinates for a point or two
 * points to construct a line. This highly depends on the elementType that is constructed. See the corresponding JXG.create*
 * methods for a list of possible parameters.
 * @param {Object} attributes An object containing the attributes to be set. This also depends on the elementType.
 * Common attributes are name, visible, strokeColor. @see GeometryElement#setProperty
 * @type Object
 * @return Reference to the created element.
 */
JXG.Board.prototype.createElement = function(elementType, parents, attributes) {
    var el, i, s;

    // CM: AW:
    if (elementType!='turtle' && (parents == null || parents.length == 0)) {  // Turtle may have no parent elements
        return null;
    }
    if (parents == null) { parents = []; }

    elementType = elementType.toLowerCase();

    if (attributes==null) {
        attributes = {};
    }
    for (i=0; i<parents.length; i++) {
        parents[i] = JXG.getReference(this, parents[i]); // TODO: should not be done for content-parameter of JXG.Text
    }

    if(JXG.JSXGraph.elements[elementType] != null) {
	if(typeof JXG.JSXGraph.elements[elementType] == 'function') {
            el = JXG.JSXGraph.elements[elementType](this, parents, attributes);
        } else {
            el = JXG.JSXGraph.elements[elementType].creator(this, parents, attributes);
        }
    } else {
        throw new Error("JSXGraph: JXG.createElement: Unknown element type given: "+elementType);
    }

    if (el==undefined) {
        //throw new Error("JSXGraph: JXG.createElement: failure creating "+elementType);
        return;
    };

    if(JXG.isArray(attributes)) {
        attributes = attributes[0];
    }

//    try {
        if(el.multipleElements) {
            for(s in el) {
                if(typeof el[s].setProperty != 'undefined')
                    el[s].setProperty(attributes);
            }
        } else {
            if(typeof el.setProperty != 'undefined')
                el.setProperty(attributes);
        }

//    } catch (e) { alert("Error setting Property:" + e); };

//    if(!JXG.isArray(el)) {  // Default way of setting attributes: strings, arrays and objects are possible
//        el.setProperty(attributes);
//    }
/* AW: Doch erstmal wieder auskommentiert
    else {                  // Setting attributes of multiple objects simultaneously.  Here, only strings are possible
        for (var s in attributes) {
            for(var i=0; i<el.length; i++) {
                if(attributes[s][i] != null) {el[i].setProperty(s+':'+attributes[s][i]);}
            }
        }
    }
*/
/*
    for (var s in attributes) {
        if(!JXG.isArray(el)) {
            el.setProperty(s+':'+attributes[s]);
        }
        else {
            for(var i=0; i<el.length; i++) {
                if(attributes[s][i] != null) {
                    el[i].setProperty(s+':'+attributes[s][i]);
                }
            }
        }
    }
*/
    this.update(el); // We start updating at the newly created element. AW
    return el;
};

/**
 * Wrapper for {@link #createElement()}.
 */
JXG.Board.prototype.create = JXG.Board.prototype.createElement;

/**
 * Delete the elements drawn as part of a trace of an element.
 */
JXG.Board.prototype.clearTraces = function() {
    var el;

    for(el in this.objects) {
        if (this.objects[el].traced)
            this.objects[el].clearTrace();
    }
    return this;
};

/**
 * Method called before a board is initialized or load from a file. Currently unused.
 * @private
 */
JXG.Board.prototype.beforeLoad = function() {
/*    if(document.getElementsByTagName("body").length > 0) {
        var divNode = document.createElement("div");
        divNode.setAttribute("id", "JXGPreLoadAnimation");
        var imgNode = document.createElement("img");
        imgNode.setAttribute("src", "./css/load.gif");
        divNode.appendChild(imgNode);
        divNode.setStyle({
                    zIndex: 999,
                    position: 'absolute',
                    left: parseInt(JXG.getStyle(this.containerObj,"left")) + (this.canvasWidth - 100)/2,
                    top: parseInt(JXG.getStyle(this.containerObj,"top")) + (this.canvasHeight - 100)/2
                });

        document.getElementsByTagName("body")[0].appendChild(divNode);
    }*/
};

/**
 * Method called after a board got initialized or load from a file. Currently unused.
 * @private
 */
JXG.Board.prototype.afterLoad = function() {
  /*  if(document.getElementsByTagName("body").length > 0) {
        document.getElementsByTagName("body")[0].removeChild(document.getElementById("JXGPreLoadAnimation"));
    }*/
};

/**
 * Stop updates of the board.
 */
JXG.Board.prototype.suspendUpdate = function() {
    this.isSuspendedUpdate = true;
};

/**
 * Enable updates of the board again.
 */
JXG.Board.prototype.unsuspendUpdate = function() {
    this.isSuspendedUpdate = false;
    this.update();
};

/**
 * Set the bounding box of the board.
 * @param {Array} New bounding box [x1,y1,x2,y2]
 * @param {Bool} keepaspectratio: optional flag
 */
JXG.Board.prototype.setBoundingBox = function(bbox,keepaspectratio) {
    if (!JXG.isArray(bbox)) return;
    var h,w,oX,oY;
    w = this.canvasWidth;
    h = this.canvasHeight;
    if (keepaspectratio) {
        this.unitX = w/(bbox[2]-bbox[0]);
        this.unitY = h/(-bbox[3]+bbox[1]);
        if (this.unitX<this.unitY) {
            this.unitY = this.unitX;
        } else {
            this.unitX = this.unitY;
        }
    } else {
        this.unitX = w/(bbox[2]-bbox[0]);
        this.unitY = h/(-bbox[3]+bbox[1]);
    }
    oX = -this.unitX*bbox[0]*this.zoomX;
    oY = this.unitY*bbox[1]*this.zoomY;
    this.origin = new JXG.Coords(JXG.COORDS_BY_SCREEN, [oX, oY], this);
    this.stretchX = this.zoomX*this.unitX;
    this.stretchY = this.zoomY*this.unitY;

    this.moveOrigin();
    return this;
};

/**
 * General purpose animation function, currently only supporting moving points from one place to another. Is faster than
 * managing the animation per point, especially if there is more than one animated point at the same time.
 */
JXG.Board.prototype.animate = function() {
    var count = 0,
        el, o, newCoords, r, p, c, 
        obj=null;

    //this.suspendUpdate();
    for(el in this.animationObjects) {
        if(this.animationObjects[el] == null)
            continue;

        count++;
        o = this.animationObjects[el];
        if(o.animationPath) {
            newCoords = o.animationPath.pop();
            if(typeof newCoords  == 'undefined') {
                delete(o.animationPath);
            } else {
                //o.setPositionByTransform(JXG.COORDS_BY_USER, newCoords[0] - o.coords.usrCoords[1], newCoords[1] - o.coords.usrCoords[2]);
                o.setPositionDirectly(JXG.COORDS_BY_USER, newCoords[0], newCoords[1]);
                //this.update(o);  // May slow down the animation, but is important 
                                 // for dependent glider objects (see tangram.html).
                                 // Otherwise the intended projection may be incorrect.
                o.prepareUpdate().update().updateRenderer();
                obj = o;
            }
        }
        if(o.animationData) {
            c = 0;
            for(r in o.animationData) {
                p = o.animationData[r].pop();
                if(typeof p == 'undefined') {
                    delete(o.animationData[p]);
                } else {
                    c++;
                    o.setProperty(r + ':' + p);
                }
            }
            if(c==0)
                delete(o.animationData);
        }

        if(typeof o.animationData == 'undefined' && typeof o.animationPath == 'undefined') {
            this.animationObjects[el] = null;
            delete(this.animationObjects[el]);
        }
    }
    //this.unsuspendUpdate();

    if(count == 0) {
        window.clearInterval(this.animationIntervalCode);
        delete(this.animationIntervalCode);
    } else {
        this.update(obj);
//	window.setTimeout('JXG.JSXGraph.boards[\'' + this.id + '\'].animate();', 35);
    }
};

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * Options object.
 * @class These are the default options of the board and
 * of all geometry elements.
 * @constructor
 */
JXG.Options = {
    /* Options that are used directly within the board class */
    fontSize : 12,
    showCopyright : true,
    showNavigation : true,
    takeSizeFromFile : false, // If true, the construction - when read from a file or string - the size of the div can be changed.
    renderer: 'svg',

    /* grid options */
    grid : {
        /* grid styles */
        hasGrid : false,
        gridX : 2,
        gridY : 2,
        gridColor : '#C0C0C0',
        gridOpacity : '0.5',
        gridDash : true,
        /* snap to grid options */
        snapToGrid : false,
        snapSizeX : 2,
        snapSizeY : 2
    },
    /* zoom options */
    zoom : {
        factor : 1.25
    },

    /* geometry element options */
    elements : {
        /* color options */
        strokeColor: '#0000ff',
        highlightStrokeColor: '#C3D9FF',
        fillColor: 'none',
        highlightFillColor: 'none',

        strokeOpacity: 1,
        highlightStrokeOpacity: 1,
        fillOpacity: 1,
        highlightFillOpacity: 1,
        strokeWidth: '2px',
	    withLabel: false,

        /*draft options */
        draft : {
            draft : false,
            color : '#565656',
            opacity : 0.8,
            strokeWidth : '1px'
        }
    },

    /* special point options */
    point : {
    	withLabel: true,
        style : 5, //1;
        fillColor : '#ff0000',
        highlightFillColor : '#EEEEEE',
        strokeWidth: '2px',
        strokeColor : '#ff0000', //'#0000ff',
        highlightStrokeColor : '#C3D9FF',
        zoom: false             // Change the point size on zoom
    },

    /* special line options */
    line : {
        firstArrow : false,
        lastArrow : false,
        straightFirst : true,
        straightLast : true,
        fillColor : '#000000',
        highlightFillColor : 'none',
        strokeColor : '#0000ff',
        highlightStrokeColor : '#888888',
        /* line ticks options */
        ticks : {
            drawLabels : true,
            drawZero : false,
            insertTicks : false,
            minTicksDistance : 50,
            maxTicksDistance : 300,
            minorHeight : 4,
            majorHeight : 10,
            minorTicks : 4,
            defaultDistance : 1
        }
    },

    /* special axis options */
    axis : {
        strokeColor : '#666666',
        highlightStrokeColor : '#888888'
    },
    
    /*special circle options */
    circle : {
        fillColor : 'none',
        highlightFillColor : 'none',
        strokeColor : '#0000ff',
        highlightStrokeColor : '#C3D9FF'
    },

    /* special conic options */
    conic : {
        fillColor : 'none',
        highlightFillColor : 'none',
        strokeColor : '#0000ff',
        highlightStrokeColor : '#C3D9FF'
    },

    /* special angle options */
    angle : {
	    withLabel:true,
        radius : 1.0,
        fillColor : '#FF7F00',
        highlightFillColor : '#FF7F00',
        strokeColor : '#FF7F00',
        fillOpacity : 0.3,
        highlightFillOpacity : 0.3
    },

    /* special arc options */
    arc : {
        firstArrow : false,
        lastArrow : false,
        fillColor : 'none',
        highlightFillColor : 'none',
        strokeColor : '#0000ff',
        highlightStrokeColor : '#C3D9FF'
    },

    /* special polygon options */
    polygon : {
        fillColor : '#00FF00',
        highlightFillColor : '#00FF00',
        fillOpacity : 0.3,
        highlightFillOpacity : 0.3
    },

    /* special sector options */
    sector : {
        fillColor: '#00FF00',
        highlightFillColor: '#00FF00',
        fillOpacity: 0.3,
        highlightFillOpacity: 0.3
    },

    /* special text options */
    text : {
        strokeColor : '#000000',
        useASCIIMathML : false,
        defaultDisplay : 'html' //'html' or 'internal'
    },

    /* special curve options */
    curve : {
        strokeWidth : '1px',
        strokeColor : '#0000ff',
        RDPsmoothing : false,    // Apply the Ramen-Douglas-Peuker algorithm
        numberPointsHigh : 1600, // Number of points on curves after mouseUp
        numberPointsLow : 400,   // Number of points on curves after mousemove
        doAdvancedPlot : true    // Use the algorithm by Gillam and Hohenwarter
                                 // It is much slower, but the result is better
    },

    /* precision options */
    precision : {
        touch    : 20,
        mouse    : 4,
        epsilon  : 0.0001,
        hasPoint : 4
    },

    // Default ordering of the layers
    layer : {
        numlayers:20, // only important in SVG
        text  : 9,
        point : 9,
        arc   : 8,
        line  : 7,
        circle: 6, 
        curve : 5,
        polygon: 4,
        sector: 3,
        angle : 2,
        grid  : 1,
        image : 0 
    }
};

/**
 * Apply the options stored in this object to all objects on the given board.
 * @param {JXG.Board} board The board to which objects the options will be applied.
 */
JXG.useStandardOptions = function(board) {
    var o = JXG.Options,
        boardHadGrid = board.hasGrid,
        el, t;

    board.hasGrid = o.grid.hasGrid;
    board.gridX = o.grid.gridX;
    board.gridY = o.grid.gridY;
    board.gridColor = o.grid.gridColor;
    board.gridOpacity = o.grid.gridOpacity;
    board.gridDash = o.grid.gridDash;
    board.snapToGrid = o.grid.snapToGrid;
    board.snapSizeX = o.grid.SnapSizeX;
    board.snapSizeY = o.grid.SnapSizeY;
    board.takeSizeFromFile = o.takeSizeFromFile;

    for(el in board.objects) {
        p = board.objects[el];
        if(p.elementClass == JXG.OBJECT_CLASS_POINT) {
            p.visProp['fillColor'] = o.point.fillColor;
            p.visProp['highlightFillColor'] = o.point.highlightFillColor;
            p.visProp['strokeColor'] = o.point.strokeColor;
            p.visProp['highlightStrokeColor'] = o.point.highlightStrokeColor;
        }
        else if(p.elementClass == JXG.OBJECT_CLASS_LINE) {
            p.visProp['fillColor'] = o.line.fillColor;
            p.visProp['highlightFillColor'] = o.line.highlightFillColor;
            p.visProp['strokeColor'] = o.line.strokeColor;
            p.visProp['highlightStrokeColor'] = o.line.highlightStrokeColor;
            for(t in p.ticks) {
                t.majorTicks = o.line.ticks.majorTicks;
                t.minTicksDistance = o.line.ticks.minTicksDistance;
                t.minorHeight = o.line.ticks.minorHeight;
                t.majorHeight = o.line.ticks.majorHeight;
            }
        }
        else if(p.elementClass == JXG.OBJECT_CLASS_CIRCLE) {
            p.visProp['fillColor'] = o.circle.fillColor;
            p.visProp['highlightFillColor'] = o.circle.highlightFillColor;
            p.visProp['strokeColor'] = o.circle.strokeColor;
            p.visProp['highlightStrokeColor'] = o.circle.highlightStrokeColor;
        }
        else if(p.type == JXG.OBJECT_TYPE_ANGLE) {
            p.visProp['fillColor'] = o.angle.fillColor;
            p.visProp['highlightFillColor'] = o.angle.highlightFillColor;
            p.visProp['strokeColor'] = o.angle.strokeColor;
        }
        else if(p.type == JXG.OBJECT_TYPE_ARC) {
            p.visProp['fillColor'] = o.arc.fillColor;
            p.visProp['highlightFillColor'] = o.arc.highlightFillColor;
            p.visProp['strokeColor'] = o.arc.strokeColor;
            p.visProp['highlightStrokeColor'] = o.arc.highlightStrokeColor;
        }
        else if(p.type == JXG.OBJECT_TYPE_POLYGON) {
            p.visProp['fillColor'] = o.polygon.fillColor;
            p.visProp['highlightFillColor'] = o.polygon.highlightFillColor;
            p.visProp['fillOpacity'] = o.polygon.fillOpacity;
            p.visProp['highlightFillOpacity'] = o.polygon.highlightFillOpacity;
        }
        else if(p.type == JXG.OBJECT_TYPE_CURVE) {
            p.visProp['strokeColor'] = o.curve.strokeColor;
        }
    }
    for(el in board.objects) {
        p = board.objects[el];
        if(p.type == JXG.OBJECT_TYPE_SECTOR) {
            p.arc.visProp['fillColor'] = o.sector.fillColor;
            p.arc.visProp['highlightFillColor'] = o.sector.highlightFillColor;
            p.arc.visProp['fillOpacity'] = o.sector.fillOpacity;
            p.arc.visProp['highlightFillOpacity'] = o.sector.highlightFillOpacity;
        }
    }

    board.fullUpdate();
    if(boardHadGrid && board.hasGrid) {
        board.renderer.removeGrid(board);
        board.renderer.drawGrid(board);
    } else if(boardHadGrid && !board.hasGrid) {
        board.renderer.removeGrid(board);
    } else if(!boardHadGrid && board.hasGrid) {
        board.renderer.drawGrid(board);
    }
};

/**
 * Converts all color values to greyscale and calls useStandardOption to put them onto the board.
 * @param {JXG.Board} board The board to which objects the options will be applied.
 * @see #useStandardOptions
 */
JXG.useBlackWhiteOptions = function(board) {
    o = JXG.Options;
    o.point.fillColor = JXG.rgb2bw(o.point.fillColor);
    o.point.highlightFillColor = JXG.rgb2bw(o.point.highlightFillColor);
    o.point.strokeColor = JXG.rgb2bw(o.point.strokeColor);
    o.point.highlightStrokeColor = JXG.rgb2bw(o.point.highlightStrokeColor);

    o.line.fillColor = JXG.rgb2bw(o.line.fillColor);
    o.line.highlightFillColor = JXG.rgb2bw(o.line.highlightFillColor);
    o.line.strokeColor = JXG.rgb2bw(o.line.strokeColor);
    o.line.highlightStrokeColor = JXG.rgb2bw(o.line.highlightStrokeColor);

    o.circle.fillColor = JXG.rgb2bw(o.circle.fillColor);
    o.circle.highlightFillColor = JXG.rgb2bw(o.circle.highlightFillColor);
    o.circle.strokeColor = JXG.rgb2bw(o.circle.strokeColor);
    o.circle.highlightStrokeColor = JXG.rgb2bw(o.circle.highlightStrokeColor);

    o.arc.fillColor = JXG.rgb2bw(o.arc.fillColor);
    o.arc.highlightFillColor = JXG.rgb2bw(o.arc.highlightFillColor);
    o.arc.strokeColor = JXG.rgb2bw(o.arc.strokeColor);
    o.arc.highlightStrokeColor = JXG.rgb2bw(o.arc.highlightStrokeColor);

    o.polygon.fillColor = JXG.rgb2bw(o.polygon.fillColor);
    o.polygon.highlightFillColor  = JXG.rgb2bw(o.polygon.highlightFillColor);

    o.sector.fillColor = JXG.rgb2bw(o.sector.fillColor);
    o.sector.highlightFillColor  = JXG.rgb2bw(o.sector.highlightFillColor);

    o.curve.strokeColor = JXG.rgb2bw(o.curve.strokeColor);
    o.grid.gridColor = JXG.rgb2bw(o.grid.gridColor);

    JXG.useStandardOptions(board);
};

/**
 * Decolorizes the given color.
 * @param {String} color HTML string containing the HTML color code.
 * @type String
 * @return Returns a HTML color string
 */
JXG.rgb2bw = function(color) {
    if(color == 'none') {
        return color;
    }
    var x, HexChars="0123456789ABCDEF", tmp, arr;
    arr = JXG.rgbParser(color);
    x = 0.3*arr[0] + 0.59*arr[1] + 0.11*arr[2];
    tmp = HexChars.charAt((x>>4)&0xf)+HexChars.charAt(x&0xf);
    color = "#" + tmp + "" + tmp + "" + tmp;
    return color;
};

/**
 * Converts the colors of the elements to how a color blind person would approximately see it. Possible
 * options are <i>protanopia</i>, <i>deuteranopia</i>, and <i>tritanopia</i>.
 * @param {JXG.Board} board The board to which objects the options will be applied.
 * @param {string} deficiency The type of deficiency which will be simulated.
 * @see #useStandardOptions
 */
JXG.simulateColorBlindness = function(board, deficiency) {
    o = JXG.Options;
    o.point.fillColor = JXG.rgb2cb(o.point.fillColor, deficiency);
    o.point.highlightFillColor = JXG.rgb2cb(o.point.highlightFillColor, deficiency);
    o.point.strokeColor = JXG.rgb2cb(o.point.strokeColor, deficiency);
    o.point.highlightStrokeColor = JXG.rgb2cb(o.point.highlightStrokeColor, deficiency);

    o.line.fillColor = JXG.rgb2cb(o.line.fillColor, deficiency);
    o.line.highlightFillColor = JXG.rgb2cb(o.line.highlightFillColor, deficiency);
    o.line.strokeColor = JXG.rgb2cb(o.line.strokeColor, deficiency);
    o.line.highlightStrokeColor = JXG.rgb2cb(o.line.highlightStrokeColor, deficiency);

    o.circle.fillColor = JXG.rgb2cb(o.circle.fillColor, deficiency);
    o.circle.highlightFillColor = JXG.rgb2cb(o.circle.highlightFillColor, deficiency);
    o.circle.strokeColor = JXG.rgb2cb(o.circle.strokeColor, deficiency);
    o.circle.highlightStrokeColor = JXG.rgb2cb(o.circle.highlightStrokeColor, deficiency);

    o.arc.fillColor = JXG.rgb2cb(o.arc.fillColor, deficiency);
    o.arc.highlightFillColor = JXG.rgb2cb(o.arc.highlightFillColor, deficiency);
    o.arc.strokeColor = JXG.rgb2cb(o.arc.strokeColor, deficiency);
    o.arc.highlightStrokeColor = JXG.rgb2cb(o.arc.highlightStrokeColor, deficiency);

    o.polygon.fillColor = JXG.rgb2cb(o.polygon.fillColor, deficiency);
    o.polygon.highlightFillColor  = JXG.rgb2cb(o.polygon.highlightFillColor, deficiency);

    o.sector.fillColor = JXG.rgb2cb(o.sector.fillColor, deficiency);
    o.sector.highlightFillColor  = JXG.rgb2cb(o.sector.highlightFillColor, deficiency);

    o.curve.strokeColor = JXG.rgb2cb(o.curve.strokeColor, deficiency);
    o.grid.gridColor = JXG.rgb2cb(o.grid.gridColor, deficiency);

    JXG.useStandardOptions(board);
};

/**
 * Decolorizes the given color.
 * @param {String} color HTML string containing the HTML color code.
 * @param {String} deficiency The type of color blindness. Possible
 * options are <i>protanopia</i>, <i>deuteranopia</i>, and <i>tritanopia</i>.
 * @type String
 * @return Returns a HTML color string
 */
JXG.rgb2cb = function(color, deficiency) {
    if(color == 'none') {
        return color;
    }

    var rgb, l, m, s, lms, tmp,
        a1, b1, c1, a2, b2, c2;
//        anchor = new Array(12), anchor_e = new Array(3);
/*
 has been required to calculate the constants for a1, ..., c2, and inflection.
*/
/* old stuff. just here for debugging purposes
    anchor[0] = 0.08008;  anchor[1]  = 0.1579;    anchor[2]  = 0.5897;
    anchor[3] = 0.1284;   anchor[4]  = 0.2237;    anchor[5]  = 0.3636;
    anchor[6] = 0.9856;   anchor[7]  = 0.7325;    anchor[8]  = 0.001079;
    anchor[9] = 0.0914;   anchor[10] = 0.007009;  anchor[11] = 0.0;

    anchor_e[0] = 0.14597772;
    anchor_e[1] = 0.12188395;
    anchor_e[2] = 0.08413913;


    document.getElementById('debug').innerHTML += 'color: ' + color;

//    document.getElementById('debug').innerHTML += 'deuteranopia<br/><br/>';    
      // find a,b,c for lam=575nm and lam=475 
      a1 = anchor_e[1] * anchor[8] - anchor_e[2] * anchor[7];
      b1 = anchor_e[2] * anchor[6] - anchor_e[0] * anchor[8];
      c1 = anchor_e[0] * anchor[7] - anchor_e[1] * anchor[6];
      a2 = anchor_e[1] * anchor[2] - anchor_e[2] * anchor[1];
      b2 = anchor_e[2] * anchor[0] - anchor_e[0] * anchor[2];
      c2 = anchor_e[0] * anchor[1] - anchor_e[1] * anchor[0];
      inflection = (anchor_e[2] / anchor_e[0]);

//    document.getElementById('debug').innerHTML += 'a1 = ' + a1 + '<br/>' + 'b1 = ' + b1 + '<br/>' + 'c1 = ' + c1 + '<br/>' + 'a2 = ' + a2 + '<br/>' + 'b2 = ' + b2 + '<br/>' + 'c2 = ' + c2 + '<br/>' + 'inflection = ' + inflection + '<br/><br/>protanopia<br/><br/>';
      // find a,b,c for lam=575nm and lam=475 
      a1 = anchor_e[1] * anchor[8] - anchor_e[2] * anchor[7];
      b1 = anchor_e[2] * anchor[6] - anchor_e[0] * anchor[8];
      c1 = anchor_e[0] * anchor[7] - anchor_e[1] * anchor[6];
      a2 = anchor_e[1] * anchor[2] - anchor_e[2] * anchor[1];
      b2 = anchor_e[2] * anchor[0] - anchor_e[0] * anchor[2];
      c2 = anchor_e[0] * anchor[1] - anchor_e[1] * anchor[0];
      inflection = (anchor_e[2] / anchor_e[1]);

//    document.getElementById('debug').innerHTML += 'a1 = ' + a1 + '<br/>' + 'b1 = ' + b1 + '<br/>' + 'c1 = ' + c1 + '<br/>' + 'a2 = ' + a2 + '<br/>' + 'b2 = ' + b2 + '<br/>' + 'c2 = ' + c2 + '<br/>' + 'inflection = ' + inflection + '<br/><br/>tritanopia<br/><br/>';
      // Set 1: regions where lambda_a=575, set 2: lambda_a=475
      a1 = anchor_e[1] * anchor[11] - anchor_e[2] * anchor[10];
      b1 = anchor_e[2] * anchor[9]  - anchor_e[0] * anchor[11];
      c1 = anchor_e[0] * anchor[10] - anchor_e[1] * anchor[9];
      a2 = anchor_e[1] * anchor[5]  - anchor_e[2] * anchor[4];
      b2 = anchor_e[2] * anchor[3]  - anchor_e[0] * anchor[5];
      c2 = anchor_e[0] * anchor[4]  - anchor_e[1] * anchor[3];
      inflection = (anchor_e[1] / anchor_e[0]);


//    document.getElementById('debug').innerHTML += 'a1 = ' + a1 + '<br/>' + 'b1 = ' + b1 + '<br/>' + 'c1 = ' + c1 + '<br/>' + 'a2 = ' + a2 + '<br/>' + 'b2 = ' + b2 + '<br/>' + 'c2 = ' + c2 + '<br/>' + 'inflection = ' + inflection;    
*/
    lms = JXG.rgb2LMS(color);
    l = lms.l; m = lms.m; s = lms.s;

    deficiency = deficiency.toLowerCase();

    switch(deficiency) {
        case "protanopia":
            a1 = -0.06150039994295001;
            b1 = 0.08277001656812001;
            c1 = -0.013200141220000003;
            a2 = 0.05858939668799999;
            b2 = -0.07934519995360001;
            c2 = 0.013289415272000003;
            inflection = 0.6903216543277437;

            tmp = s/m;
            if (tmp < inflection)
                l = -(b1 * m + c1 * s) / a1;
            else
                l = -(b2 * m + c2 * s) / a2;
            break;
        case "tritanopia":
            a1 = -0.00058973116217;
            b1 = 0.007690316482;
            c1 = -0.01011703519052;
            a2 = 0.025495080838999994;
            b2 = -0.0422740347;
            c2 = 0.017005316784;
            inflection = 0.8349489908460004;

            tmp = m / l;
            if (tmp < inflection)
              s = -(a1 * l + b1 * m) / c1;
            else
              s = -(a2 * l + b2 * m) / c2;
            break;
        default:
            a1 = -0.06150039994295001;
            b1 = 0.08277001656812001;
            c1 = -0.013200141220000003;
            a2 = 0.05858939668799999;
            b2 = -0.07934519995360001;
            c2 = 0.013289415272000003;
            inflection = 0.5763833686400911;

            tmp = s/l;
            if(tmp < inflection)
                m = -(a1 * l + c1 * s) / b1;
            else
                m = -(a2 * l + c2 * s) / b2;
            break;
    }

    rgb = JXG.LMS2rgb(l, m, s);

    var HexChars="0123456789ABCDEF";
    tmp = HexChars.charAt((rgb.r>>4)&0xf)+HexChars.charAt(rgb.r&0xf);
    color = "#" + tmp;
    tmp = HexChars.charAt((rgb.g>>4)&0xf)+HexChars.charAt(rgb.g&0xf);
    color += tmp;
    tmp = HexChars.charAt((rgb.b>>4)&0xf)+HexChars.charAt(rgb.b&0xf);
    color += tmp;

    return color;
};

/**
 * Load options from a file using FileReader
 * @param fileurl {String} URL to .json-file containing style information
 * @param apply {bool} <tt>true</tt> when options in file should be applied to board after being loaded.
 * @param board {JXG.Board} The board the options should be applied to.
 */
JXG.loadOptionsFromFile = function(fileurl, applyTo, board) {
   this.cbp = function(t) {
      this.parseString(t, applyTo, board);
   };
   this.cb = JXG.bind(this.cbp,this);

   JXG.FileReader.parseFileContent(fileurl, this.cb, 'raw');
};

/**
 * Apply options given as a string to a board.
 * @param text {String} Options given as a string in .json-Format
 * @param apply {bool} <tt>true</tt> if the options should be applied to all objects on the board.
 * @param board {JXG.Board} The board the options should be applied to.
 */
JXG.parseOptionsString = function(text, applyTo, board) {
   var newOptions = '';

   if(text != '') {
      newOptions = eval("(" + text + ")");
   }
   else
      return;

   var maxDepth = 10;
   var applyOption = function (base, option, depth) {
      if(depth==10)
         return;
      depth++;

      for(var key in option) {
         if((JXG.isNumber(option[key])) || (JXG.isArray(option[key])) || (JXG.isString(option[key])) || (option[key]==true) || (option[key]==false)) {
            base[key] = option[key];
         }
         else {
            applyOption(base[key], option[key], depth);
         }
      }
   };

   applyOption(this, newOptions, 0);

   if(applyTo && typeof board != 'undefined') {
       JXG.useStandardOptions(board);
   }
};

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview Class Geonext is defined in this file. Geonext controls all boards.
 * It has methods to create, save, load and free boards.
 * @author graphjs
 * @version 0.81
 */

/**
 * Constructs a new Geonext singleton object.
 * @class This is the Geonext class. It stores all properties required
 * to load, save, create and free a board.
 * @constructor
 * @param {String} forceRenderer If a specific renderer should be chosen. Possible values are 'vml', 'svg', 'silverlight'
 */
JXG.JSXGraph = new function () {
    var ie, opera, i, arr;
    this.licenseText = 'JSXGraph v0.81rc1 Copyright (C) see http://jsxgraph.org';

    /**
            * Stores the renderer that is used to draw the board.
            * @type String
            */
    this.rendererType = '';

    /**
            * Associative array that keeps all boards.
            * @type Object
            */
    this.boards = {};

    /**
            * Associative array that keeps all registered geometry elements
            * @type Object
            */
    this.elements = {};

    if( (typeof forceRenderer == 'undefined') || (forceRenderer == null) || (forceRenderer == '') ) {
        /* Determine the users browser */
        ie = navigator.appVersion.match(/MSIE (\d\.\d)/);
        opera = (navigator.userAgent.toLowerCase().indexOf("opera") != -1);

        /* and set the rendererType according to the browser */
        if ((!ie) || (opera)) {
            //this.rendererType = 'svg';
            JXG.Options.renderer = 'svg';
        }
        else {
            //if(Silverlight.available)
            //    this.rendererType = 'silverlight';
            //else
                //this.rendererType = 'vml';
                JXG.Options.renderer = 'vml';
                function MouseMove(e) { //Magic!
                  document.body.scrollLeft;
                  document.body.scrollTop;
                }
                document.onmousemove = MouseMove;
        }
    } else {
        /* the user has chosen a specific renderer */
        this.rendererType = forceRenderer;
    }

    /* Load the source files for the renderer */
    //JXG.rendererFiles[this.rendererType].split(',').each( function(include) { JXG.require(JXG.requirePath+include+'.js'); } );
    arr = JXG.rendererFiles[JXG.Options.renderer].split(',');
    for (i=0;i<arr.length;i++) ( function(include) { JXG.require(JXG.requirePath+include+'.js'); } )(arr[i]);


    /**
            * Initialise a new board.
            * @param {String} box Html-ID to the Html-element in which the board is painted.
            * @return {JXG.Board} Reference to the created board.
            */
    this.initBoard = function (box, attributes) {
        // Create a new renderer
        var renderer,
            originX, originY, unitX, unitY,
            w, h, dimensions,
            bbox,
            zoomfactor, zoomX, zoomY,
            showCopyright, showNavi,
            board;

        dimensions = JXG.getDimensions(box);
        if (typeof attributes == 'undefined') {
            attributes = {};
        }
        if (typeof attributes["boundingbox"] != 'undefined') {
            bbox = attributes["boundingbox"];
            w = parseInt(dimensions.width);
            h = parseInt(dimensions.height);
            if (attributes["keepaspectratio"]) {
            /**
                                * If the boundingbox attribute is given and the ratio of height and width of the sides defined by the bounding box and
                                * the ratio of the dimensions of the div tag which contains the board do not coincide,
                                * then the smaller side is chosen.
                                */
                unitX = w/(bbox[2]-bbox[0]);
                unitY = h/(-bbox[3]+bbox[1]);
                if (unitX<unitY) {
                    unitY = unitX;
                } else {
                    unitX = unitY;
                }
            } else {
                unitX = w/(bbox[2]-bbox[0]);
                unitY = h/(-bbox[3]+bbox[1]);
            }
            originX = -unitX*bbox[0];
            originY = unitY*bbox[1];
        } else {
            originX = ( (typeof attributes["originX"]) == 'undefined' ? 150 : attributes["originX"]);
            originY = ( (typeof attributes["originY"]) == 'undefined' ? 150 : attributes["originY"]);
            unitX = ( (typeof attributes["unitX"]) == 'undefined' ? 50 : attributes["unitX"]);
            unitY = ( (typeof attributes["unitY"]) == 'undefined' ? 50 : attributes["unitY"]);
        }
        zoomfactor = ( (typeof attributes["zoom"]) == 'undefined' ? 1.0 : attributes["zoom"]);
        zoomX = zoomfactor*( (typeof attributes["zoomX"]) == 'undefined' ? 1.0 : attributes["zoomX"]);
        zoomY = zoomfactor*( (typeof attributes["zoomY"]) == 'undefined' ? 1.0 : attributes["zoomY"]);

        // ??? if (typeof attributes["showcopyright"] != 'undefined') attributes["showCopyright"] = attributes["showcopyright"];
        showCopyright = ( (typeof attributes["showCopyright"]) == 'undefined' ? JXG.Options.showCopyright : attributes["showCopyright"]);

        if(JXG.Options.renderer == 'svg') {
            renderer = new JXG.SVGRenderer(document.getElementById(box));
        } else if(JXG.Options.renderer == 'vml') {
            renderer = new JXG.VMLRenderer(document.getElementById(box));
        } else {
            renderer = new JXG.SilverlightRenderer(document.getElementById(box), dimensions.width, dimensions.height);
        }

        board = new JXG.Board(box, renderer, '', [originX, originY], 1.0, 1.0, unitX, unitY, dimensions.width, dimensions.height,showCopyright);
        this.boards[board.id] = board;
        // board.initGeonextBoard();  // Contsruct "Ursprung" and other elements.
        board.initInfobox();

        if((typeof attributes["axis"] != 'undefined') && attributes["axis"]) {
        	board.defaultAxes = {};
            board.defaultAxes.x = board.create('axis', [[0,0], [1,0]], {});
            board.defaultAxes.y = board.create('axis', [[0,0], [0,1]], {});
        }

        if ((typeof attributes["grid"] != 'undefined') && attributes["grid"]) {
            board.renderer.drawGrid(board);
        }

        if (typeof attributes["shownavigation"] != 'undefined') attributes["showNavigation"] = attributes["shownavigation"];
        showNavi = ( (typeof attributes["showNavigation"]) == 'undefined' ? board.options.showNavigation : attributes["showNavigation"]);
        if (showNavi) {
            board.renderer.drawZoomBar(board);
        }

        return board;
    };

    /**
     * Load a board from a file of format GEONExT or Intergeo.
     * @param {String} box Html-ID to the Html-element in which the board is painted.
     * @param {String} file Url to the geonext-file.
     * @param {String} string containing the file format: 'Geonext' or 'Intergeo'.
     * @return {JXG.Board} Reference to the created board.
     * @see JXG.GeonextReader
     */
    this.loadBoardFromFile = function (box, file, format) {
        var renderer, board, dimensions;

        if(JXG.Options.renderer == 'svg') {
            renderer = new JXG.SVGRenderer(document.getElementById(box));
        } else {
            renderer = new JXG.VMLRenderer(document.getElementById(box));
        }
        //var dimensions = document.getElementById(box).getDimensions();
        dimensions = JXG.getDimensions(box);

        /* User default parameters, in parse* the values in the gxt files are submitted to board */
        board = new JXG.Board(box, renderer, '', [150, 150], 1.0, 1.0, 50, 50, dimensions.width, dimensions.height);
        board.initInfobox();
        board.beforeLoad();
        JXG.FileReader.parseFileContent(file, board, format);
        if(board.options.showNavigation) {
            board.renderer.drawZoomBar(board);
        }
        this.boards[board.id] = board;
        return board;
    };

    /**
     * Load a board from a base64 encoded string containing a GEONExT or Intergeo construction.
     * @param {String} box Html-ID to the Html-element in which the board is painted.
     * @param {String} string base64 encoded string.
     * @param {String} string containing the file format: 'Geonext' or 'Intergeo'.
     * @return {JXG.Board} Reference to the created board.
     * @see JXG.GeonextReader
     */
    this.loadBoardFromString = function(box, string, format) {
        var renderer, dimensions, board;

        if(JXG.Options.renderer == 'svg') {
            renderer = new JXG.SVGRenderer(document.getElementById(box));
        } else {
            renderer = new JXG.VMLRenderer(document.getElementById(box));
        }
        //var dimensions = document.getElementById(box).getDimensions();
        dimensions = JXG.getDimensions(box);

        /* User default parameters, in parse* the values in the gxt files are submitted to board */
        board = new JXG.Board(box, renderer, '', [150, 150], 1.0, 1.0, 50, 50, dimensions.width, dimensions.height);
        board.initInfobox();
        board.beforeLoad();

        JXG.FileReader.parseString(string, board, format, true);
        if (board.options.showNavigation) {
            board.renderer.drawZoomBar(board);
        }

        this.boards[board.id] = board;
        return board;
    };

    /**
     * Free a board.
     * @param {String} box Html-ID to the Html-element in which the board was painted.
     */
    this.freeBoard = function (board) {
        var el;

        if(typeof(board) == 'string') {
            board = this.boards[board];
        }

        // Remove the event listeners
        JXG.removeEvent(document, 'mousedown', board.mouseDownListener, board);
        JXG.removeEvent(document, 'mouseup', board.mouseUpListener,board);
        JXG.removeEvent(board.containerObj, 'mousemove', board.mouseMoveListener, board);

        // Remove all objects from the board.
        for(el in board.objects) {
            board.removeObject(board.objects[el]);
        }

        // Remove all the other things, left on the board
        board.containerObj.innerHTML = '';

        // Tell the browser the objects aren't needed anymore
        for(el in board.objects) {
            delete(board.objects[el]);
        }

        // Free the renderer and the algebra object
        delete(board.renderer);
        delete(board.algebra);

        // Finally remove the board itself from the boards array
        delete(this.boards[board.id]);
    };

    this.registerElement = function (element, creator) {
        element = element.toLowerCase();
        this.elements[element] = creator;

        if(JXG.Board.prototype['_' + element])
        	throw new Error("JSXGraph: Can't create wrapper method in JXG.Board because member '_" + element + "' already exists'");
        JXG.Board.prototype['_' + element] = function (parents, attributes) {
        	return this.create(element, parents, attributes);
        };

    };

    this.unregisterElement = function (element) {
        delete (this.elements[element.toLowerCase()]);
        delete (JXG.Board.prototype['_' + element.toLowerCase()]);
    };
};

/**
 * Parameter magic: object may be a string containing the name or id of the object or
 * even the object itself, this function gets a returns to the object. Order: id/object, name.
 * @param {JXG.Board} board Reference to the board the object belongs to.
 * @param {String,Object} object String or reference to the object the reference is needed.
 * @return {Object} Reference to the object given in parameter object
 */
JXG.getReference = function(board, object) {
    if(typeof(object) == 'string') {
        if(board.objects[object] != null) { // Search by ID
            object = board.objects[object];
        } else if (board.elementsByName[object] != null) { // Search by name
            object = board.elementsByName[object];
        }
    }

    return object;
};

JXG.isString = function(obj) {
    return typeof obj == "string";
};

JXG.isNumber = function(obj) {
    return typeof obj == "number";
};

JXG.isFunction = function(obj) {
    return typeof obj == "function";
};

JXG.isArray = function(obj) {
    // Borrowed from prototype.js
    return obj != null && typeof obj == "object" && 'splice' in obj && 'join' in obj;
};

JXG.isPoint = function(p) {
    if(typeof p == 'object') {
        return (p.elementClass == JXG.OBJECT_CLASS_POINT);
    }

    return false;
};

/**
 * Converts a string containing either <strong>true</strong> or <strong>false</strong> into a boolean value.
 * @param s String containing either <strong>true</strong> or <strong>false</strong>.
 * @return String typed boolean value converted to boolean.
 */
JXG.str2Bool = function(/** string */ s) /** boolean */ {
    if (s==undefined || s==null) {
        return true;
    }
    if (typeof s == 'boolean') { 
        return s;
    }
    if (s.toLowerCase()!='true') {
        return false;
    } else {
        return true;
    }
};

JXG._board = function(box, attributes) {
	return JXG.JSXGraph.initBoard(box, attributes);
};

/**
  * Convert String, number or function to function.
  * This method is used in Transformation.js
  */
JXG.createEvalFunction = function(board,param,n) {
    // convert GEONExT syntax into function
    var f = [], i, str;

    for (i=0;i<n;i++) {
        if (typeof param[i] == 'string') {
            str = board.algebra.geonext2JS(param[i]);
            str = str.replace(/this\.board\./g,'board.');
            f[i] = new Function('','return ' + (str) + ';');
        }
    }
    return function(k) {
        var a = param[k];
        if (typeof a == 'string') {
            return f[k]();
        } else if (typeof a=='function') {
            return a();
        } else if (typeof a=='number') {
            return a;
        }
        return 0;
    };
};

/**
  * Convert String, number or function to function.
  **/
JXG.createFunction = function(term,board,variableName,evalGeonext) {
    var newTerm;

    if ((evalGeonext==null || evalGeonext==true) && JXG.isString(term)) {
        // Convert GEONExT syntax into  JavaScript syntax
        newTerm = board.algebra.geonext2JS(term);
        return new Function(variableName,'return ' + newTerm + ';');
    } else if (JXG.isFunction(term)) {
        return term;
    } else if (JXG.isNumber(term)) {
        return function() { return term; };
    } else if (JXG.isString(term)) {        // In case of string function like fontsize
        return function() { return term; };
    }
    return null;
};

/*
JXG.checkParameter = function(board, parameter, input, output) {
    var r;
    if (input=='point') {
        if (JXG.isPoint(input) && output=='point') { return parameter; }
        if (JXG.isString(input) && output=='point') { 
            r = JXG.getReference(board,parameter);
            if (JXG.isString(r)) { return false; } else { return r; }
        }
    } else if (input=='array') {
        if (JXG.isArray(input) && output=='point') { 
            return = board.create('point', parameter, {visible:false,fixed:true});
        }
    } else if (input=='line') {
...    
    }
}

JXG.readParameter = function(board, parameter, input, output) {
    var i, j, lenOut = output.length, 
        len, result;

    if (lenOut==1) {
        len = input.length;
        for (j=0;j<len;j++) {
            result = JXG.checkParameter(board, parameter, input[j], output[0]);
            if (result!=false) return result;
        }
    } else {
        for (i=0;i<lenOut;i++) {
            len = input[i].length;
            for (j=0;j<len;j++) {
                result = JXG.checkParameter(board, parameter, input[i][j], output[i]);
                if (result!=false) return result;
            }
        }
    }
    return false;
};
*/

JXG.readOption = function(options, eltype, key) {
    var val = options.elements[key];
    if (typeof options[eltype][key]!='undefined') val = options[eltype][key];
    return val;
};

JXG.checkAttributes = function(atts, keyvaluepairs) {
    var key;
    if (atts==null) { atts = {}; }
    for (key in keyvaluepairs) {
        if(atts[key] == null || typeof atts[key] == 'undefined') {
            atts[key] = keyvaluepairs[key];
        }
    }
    return atts;
};

JXG.getDimensions = function(elementId) {
    var element, display, els, originalVisibility, originalPosition,
        originalDisplay, originalWidth, originalHeight;

    // Borrowed from prototype.js
    element = document.getElementById(elementId);
    if (element==null) {
        throw new Error("\nJSXGraph: HTML container element '" + (elementId) + "' not found.");
    }

    display = element.style['display'];
    if (display != 'none' && display != null) {// Safari bug
        return {width: element.offsetWidth, height: element.offsetHeight};
    }

    // All *Width and *Height properties give 0 on elements with display none,
    // so enable the element temporarily
    els = element.style;
    originalVisibility = els.visibility;
    originalPosition = els.position;
    originalDisplay = els.display;
    els.visibility = 'hidden';
    els.position = 'absolute';
    els.display = 'block';

    originalWidth = element.clientWidth;
    originalHeight = element.clientHeight;
    els.display = originalDisplay;
    els.position = originalPosition;
    els.visibility = originalVisibility;
    return {width: originalWidth, height: originalHeight};
};

/**
  * addEvent.
  */
JXG.addEvent = function( obj, type, fn, owner ) {
    owner['x_internal'+type] = function() {return fn.apply(owner,arguments);};
    if (typeof obj.addEventListener!='undefined') { // Non-IE browser
        obj.addEventListener(type, owner['x_internal'+type], false);
    } else {  // IE
        obj.attachEvent('on'+type, owner['x_internal'+type]);
    }
};

/**
  * removeEvent.
  */
JXG.removeEvent = function( obj, type, fn, owner ) {
    try {
        if (typeof obj.addEventListener!='undefined') { // Non-IE browser
            obj.removeEventListener(type, owner['x_internal'+type], false);
        } else {  // IE
            obj.detachEvent('on'+type, owner['x_internal'+type]);
        }
    } catch(e) {
        //document.getElementById('debug').innerHTML += 'on'+type+': ' + owner['x_internal'+type]+'<br>\n';
    }
};

JXG.bind = function(fn, owner ) {
    return function() {
        return fn.apply(owner,arguments);
    };
};

/**
  * getPosition: independent from prototype and jQuery
  */
JXG.getPosition = function (Evt) {
    var posx = 0,
        posy = 0,
        Evt;

    if (!Evt) {
        Evt = window.event;
    }

    if (Evt.pageX || Evt.pageY)     {
        posx = Evt.pageX;
        posy = Evt.pageY;
    }
    else if (Evt.clientX || Evt.clientY)    {
        posx = Evt.clientX + document.body.scrollLeft + document.documentElement.scrollLeft;
        posy = Evt.clientY + document.body.scrollTop + document.documentElement.scrollTop;
    }
    return [posx,posy];
};

/**
  * getOffset: Abstraction layer for Prototype.js and jQuery
  */
JXG.getOffset = function (obj) {
    var o=obj,
        l=o.offsetLeft,
        t=o.offsetTop;

    while(o=o.offsetParent) {
        l+=o.offsetLeft;
        t+=o.offsetTop;
        if(o.offsetParent) {
            l+=o.clientLeft;
            t+=o.clientTop;
        }
    }
    return [l,t];
};

/*
JXG.getOffset = function (obj) {
    var o;

    if (typeof Prototype!='undefined' && typeof Prototype.Browser!='undefined') { // Prototype lib
        return Element.cumulativeOffset(obj);
    } else {                         // jQuery
        o = $(obj).offset();
        return [o.left,o.top];
    }
};
*/

/**
  * getStyle: Abstraction layer for Prototype.js and jQuery
  * Now independent from Prototype  and jQuery
  */
JXG.getStyle = function (obj, stylename) {
    return obj.style[stylename];
/*
    if (typeof Prototype!='undefined' && typeof Prototype.Browser!='undefined') { // Prototype lib
        return $(obj).getStyle(stylename);
    } else {
        if (typeof $(obj).attr(stylename)!='undefined') {
            return $(obj).attr(stylename);
        } else {
            return $(obj).css(stylename);
        }
    }
*/
};

JXG.keys = function(object) {
    var keys = [], property;

    for (property in object) {
        keys.push(property);
    }
    return keys;
};

JXG.escapeHTML = function(str) {
    return str.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
};

JXG.unescapeHTML = function(str) {
    return str.replace(/<\/?[^>]+>/gi, '').replace(/&amp;/g,'&').replace(/&lt;/g,'<').replace(/&gt;/g,'>');
};

/**
 * This outputs an object with a base class reference to the given object. This is useful if
 * you need a copy of an e.g. attributes object and want to overwrite some of the attributes
 * without changing the original object.
 * @param {Object} obj Object to be embedded.
 * @type Object
 * @return An object with a base class reference to <tt>obj</tt>.
 */
JXG.clone = function(obj) {
    var cObj = {};
    cObj.prototype = obj;
    return cObj;
};

/**
 * Outputs a deep copy of an existing object and not only a flat copy.
 * @param {Object} obj Object to be copied.
 * @type Object
 * @return Deep copy of given object.
 */
JXG.deepCopy = function(obj) {
    var c, i, prop, j;

    if (typeof obj !== 'object' || obj == null) {
        return obj;
    }
    if (this.isArray(obj)) {
        c = [];
        for (i=0; i<obj.length; i++) {
            prop = obj[i];
            if (typeof prop == 'object') {
                if (this.isArray(prop)) {
                    c[i] = [];
                    for (j = 0; j < prop.length; j++) {
                        if (typeof prop[j] != 'object') {
                            c[i].push(prop[j]);
                        } else {
                            c[i].push(this.deepCopy(prop[j]));
                        }
                    }
                } else {
                    c[i] = this.deepCopy(prop);
                }
            } else {
                c[i] = prop;
            }
        }
    } else {
        c = {};
        for (i in obj) {
            prop = obj[i];
            if (typeof prop == 'object') {
                if (this.isArray(prop)) {
                    c[i] = [];
                    for (j = 0; j < prop.length; j++) {
                        if (typeof prop[j] != 'object') {
                            c[i].push(prop[j]);
                        } else {
                            c[i].push(this.deepCopy(prop[j]));
                        }
                    }
                } else {
                    c[i] = this.deepCopy(prop);
                }
            } else {
                c[i] = prop;
            }
        }
    }
    return c;
};

/**
 * Embeds an existing object into another one just like {@link #clone} and copies the contents of the second object
 * to the new one. Warning: The copied properties of obj2 are just flat copies.
 * @param {Object} obj Object to be copied.
 * @param {Object} obj2 Object with data that is to be copied to the new one as well.
 * @type Object
 * @return Copy of given object including some new/overwritten data from obj2.
 */
JXG.cloneAndCopy = function(obj, obj2) {
    var cObj = {}, r;
    cObj.prototype = obj;
    for(r in obj2)
        cObj[r] = obj2[r];

    return cObj;
};

JXG.toJSON = function(obj) {
    switch (typeof obj) {
        case 'object':
            if (obj) {
                var list = [];
                if (obj instanceof Array) {
                    for (var i=0;i < obj.length;i++) {
                        list.push(JXG.toJSON(obj[i]));
                    }
                    return '[' + list.join(',') + ']';
                } else {
                    for (var prop in obj) {
                        list.push('"' + prop + '":' + JXG.toJSON(obj[prop]));
                    }
                    return '{' + list.join(',') + '}';
                }
            } else {
                return 'null';
            }
        case 'string':
            return '"' + obj.replace(/(["'])/g, '\\$1') + '"';
        case 'number':
        case 'boolean':
            return new String(obj);
    }
};

JXG.capitalize = function(str) {
    return str.charAt(0).toUpperCase() + str.substring(1).toLowerCase();
};

/**
 * Copyright 2009 Nicholas C. Zakas. All rights reserved.
 * MIT Licensed
 * param array items to do
 * param function function that is applied for everyl array item
 * param object context meaning of this in function process
 * param function callback function called after the last array element has been processed.
**/
JXG.timedChunk = function(items, process, context, callback) {
    var todo = items.concat();   //create a clone of the original
    setTimeout(function(){
        var start = +new Date();
        do {
            process.call(context, todo.shift());
        } while (todo.length > 0 && (+new Date() - start < 300));
        if (todo.length > 0){
            setTimeout(arguments.callee, 1);
        } else {
            callback(items);
        }
    }, 1);
};

JXG.trimNumber = function(str) {
	str = str.replace(/^0+/, "");
	str = str.replace(/0+$/, "");
	if(str[str.length-1] == '.' || str[str.length-1] == ',') {
		str = str.slice(0, -1);
	}
    if(str[0] == '.' || str[0] == ',') {
        str = "0" + str;
    }
	
	return str;
};

JXG.trim = function(str) {
	str = str.replace(/^w+/, "");
	str = str.replace(/w+$/, "");
	
	return str;
};

/*
JXG.isSilverlightInstalled = function() {
    var isInstalled = false,
        activeX, tryOtherBrowsers, slPlugin;

    try
    {
        activeX = null;
        tryOtherBrowsers = false;

        if (window.ActiveXObject)
        {
            try
            {
                activeX = new ActiveXObject('AgControl.AgControl');
                isInstalled = true;
                activeX = null;
            }
            catch (e)
            {
                tryOtherBrowsers = true;
            }
        }
        else
        {
            tryOtherBrowsers = true;
        }
        if (tryOtherBrowsers)
        {
            slPlugin = navigator.plugins["Silverlight Plug-In"];
            if (slPlugin)
            {
                    isInstalled = true;
            }
        }
    }
    catch (e)
    {
        isInstalled = false;
    }

    return isInstalled;
};
*/

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/
JXG.OBJECT_TYPE_ARC  = 0x4F544143;                 // Hex fuer OTAC = Object Type ArC
JXG.OBJECT_TYPE_ARROW  = 0x4F544157;                 // Hex fuer OTAW = Object Type ArroW
JXG.OBJECT_TYPE_AXIS  = 0x4F544158;                 // Hex fuer OTAX = Object Type AXis
JXG.OBJECT_TYPE_TICKS  = 0x4F545458;                 // Hex fuer OTTX = Object Type TiX
JXG.OBJECT_TYPE_CIRCLE  = 0x4F54434C;                 // Hex fuer OTCC = Object Type CirCle
JXG.OBJECT_TYPE_CURVE  = 0x4F544750;                 // Hex fuer OTGP = Object Type GraphPlot
JXG.OBJECT_TYPE_GLIDER  = 0x4F54474C;                 // Hex fuer OTGL = Object Type GLider
JXG.OBJECT_TYPE_IMAGE  = 0x4F54524D;                 // Hex fuer OTIM = Object Type IMage
JXG.OBJECT_TYPE_LINE  = 0x4F544C4E;                 // Hex fuer OTLN = Object Type LiNe
JXG.OBJECT_TYPE_POINT  = 0x4F545054;                 // Hex fuer OTPT = Object Type PoinT
JXG.OBJECT_TYPE_SLIDER = 0x4F545344;                 // Hex fuer OTSD = Object Type SliDer
JXG.OBJECT_TYPE_CAS    = 0x4F544350;                 // Hex fuer OTCP = Object Type CasPoint
JXG.OBJECT_TYPE_POLYGON  = 0x4F545059;                 // Hex fuer OTPY = Object Type PolYgon
JXG.OBJECT_TYPE_SECTOR  = 0x4F545343;                 // Hex fuer OTSC = Object Type SeCtor
JXG.OBJECT_TYPE_TEXT  = 0x4F545445;                 // Hex fuer OTTE = Object Type TextElement
JXG.OBJECT_TYPE_ANGLE = 0x4F544147;                 // Hex fuer OTAG = Object Type AnGle
JXG.OBJECT_TYPE_INTERSECTION = 0x4F54524E;          // Hex fuer OTIN = Object Type INtersection
JXG.OBJECT_TYPE_TURTLE = 0x4F5455;                 // Hex fuer OTTU = Object Type TUrtle

JXG.OBJECT_CLASS_POINT = 1;
JXG.OBJECT_CLASS_LINE = 2;
JXG.OBJECT_CLASS_CIRCLE = 3;
JXG.OBJECT_CLASS_CURVE = 4;
JXG.OBJECT_CLASS_AREA = 5;
JXG.OBJECT_CLASS_OTHER = 6;

/**
 * Constructs a new GeometryElement object.
 * @class This is the basic class for geometry elements like points, circles and lines.
 * @constructor
 * of identical elements on the board. Is not yet implemented for all elements, only points, lines and circle can be traced.
 */
JXG.GeometryElement = function() {
    /**
     * Reference to board where the element is drawn
     * @type JXG.Board
     * @default null
     * @see JXG.Board
     * @private
     */
    this.board = null;

    /**
     * Unique identifier for the element. Equivalent to id-attribute of renderer element.
     * @type String
     * @default empty string
     * @private
     */
    this.id = '';

    /**
     * Controls if updates are necessary
     * @type bool
     * @default true
     * @private
     */
    this.needsUpdate = true;

    /**
     * Not necessarily unique name for the element.
     * @type String
     * @default Name generated by {@link JXG.Board#generateName}.
     * @see JXG.Board#generateName
     */
    this.name = '';

    /**
     * An associative array containing all visual properties.
     * @type Object
     * @default empty object
     * @private
     */
    this.visProp = {};

    JXG.clearVisPropOld(this); // create this.visPropOld and set default values

    /**
     * If element is in two dimensional real space this is true, else false.
     * @type boolean
     * @default true
     * @private
     */
    this.isReal = true;

    /**
     * Determines the elements border-style.
     * Possible values are:
     * <ul><li>0 for a solid line</li>
     * <li>1 for a dotted line</li>
     * <li>2 for a line with small dashes</li>
     * <li>3 for a line with medium dashes</li>
     * <li>4 for a line with big dashes</li>
     * <li>5 for a line with alternating medium and big dashes and large gaps</li>
     * <li>6 for a line with alternating medium and big dashes and small gaps</li></ul>
     * @type number
     * @name JXG.GeometryElement#dash
     * @default 0
     */
    this.visProp['dash'] = 0;

    /**
     * Stores all dependent objects to be updated when this point is moved.
     * @type Object
     * @private
     */
    this.childElements = {};

    /**
     * If element has a label subelement then this property will be set to true.
     * @type boolean
     * @default false
     * @private
     */
    this.hasLabel = false;

    /**
     * display layer which will conting the element.
     * Controlled in JXG.Options.
     */
    this.layer = 9;

    /**
     * Stores all Intersection Objects which in this moment are not real and
     * so hide this element.
     * @type object
     * @private
     */
    this.notExistingParents = {};

    /**
     * If true the element will be traced, i.e. on every movement the element will be copied
     * to the background. Use {@link #clearTrace} to delete the trace elements.
     * @see #clearTrace
     * @see #traces
     * @see #numTraces
     * @type boolean
     * @default false
     * @name JXG.GeometryElement#trace
     */
    this.traced = false;

    /**
     * Keeps track of all objects drawn as part of the trace of the element.
     * @see #traced
     * @see #clearTrace
     * @see #numTraces
     * @type Object
     * @private
     */
    this.traces = {};

    /**
     * Counts the number of objects drawn as part of the trace of the element.
     * @see #traced
     * @see #clearTrace
     * @see #traces
     * @type number
     * @private
     */
    this.numTraces = 0;

    /**
     * Stores the  transformations which are applied during update in an array
     * @type Array
     * @see JXG.Transformation
     * @private
     */
    this.transformations = [];

    /** TODO
     * @type TODO
     * @default null
     * @private
     */
    this.baseElement = null;

    /**
     * Elements depending on this element are stored here.
     * @type object
     * @private
     */
    this.descendants = {};

    /**
     * Elements on which this elements depends on are stored here.
     * @type object
     * @private
     */
    this.ancestors = {};

    /**
     * Stores variables for symbolic computations
     * @type Object
     * @private
     */
    this.symbolic = {};

    /**
     * [c,b0,b1,a,k,r,q0,q1]
     *
     * See
     * A.E. Middleditch, T.W. Stacey, and S.B. Tor:
     * "Intersection Algorithms for Lines and Circles",
     * ACM Transactions on Graphics, Vol. 8, 1, 1989, pp 25-40.
     *
     * The meaning of the parameters is:
     * Circle: points p=[p0,p1] on the circle fulfill
     *  a&lt;p,p&gt; + &lt;b,p&gt; + c = 0
     * For convenience we also store
     *  r: radius
     *  k: discriminant = sqrt(&lt;b,b&gt;-4ac)
     *  q=[q0,q1] center
     *
     * Points have radius = 0.
     * Lines have radius = infinity.
     * b: normalized vector, representing the direction of the line.
     *
     * Should be put into Coords, when all elements possess Coords.
     * @type array
     * @default [1, 0, 0, 0, 1, 1, 0, 0]
     * @private
     */
    this.stdform = [1,0,0,0,1, 1,0,0];

    /**
     * Quadratic form representation of circles (and conics)
     */
    this.quadraticform = [[1,0,0],[0,1,0],[0,0,1]];

    /**
     * If this is set to true, the element is updated in every update
     * call of the board. If set to false, the element is updated only after
     * zoom events or more generally, when the bounding box has been changed.
     * Examples for the latter behaviour should be axes.
     *
     * @type boolean
     * @default true
     * @private
     */
    this.needsRegularUpdate = true;

};

/**
 * Initializes board, id and name which cannot be initialized properly in the constructor.
 * @param {String,JXG.Board} board The board the new point is drawn on.
 * @param {String} id Unique identifier for the point. If null or an empty string is given,
 *  an unique id will be generated by Board
 * @param {String} name Not necessarily unique name for the point. If null or an
 *  empty string is given, an unique name will be generated
 * @private
 */
JXG.GeometryElement.prototype.init = function(board, id, name) {
    /*
     * Parameter magic, if board is a string, assume it is an if of an object of
     * type Board an get the boards reference.
     */
    if (typeof(board) == 'string') {
        board = JXG.JSXGraph.boards[board];
    }

    /* already documented in constructor */
    this.board = board;

    /* already documented in constructor */
    this.id = id;

    /* If name is not set or null or even undefined, generate an unique name for this object */
    if ( /*(name != '') &&*/ (name != null) && (typeof name != 'undefined') ) {
        name = name;
    } else {
        name = this.board.generateName(this);
    }
    this.board.elementsByName[name] = this;

    /* already documented in constructor */
    this.name = name;

    /**
     * The stroke color of the given geometry element.
     * @type string
     * @name JXG.GeometryElement#strokeColor
     * @see #highlightStrokeColor
     * @see #strokeWidth
     * @see #strokeOpacity
     * @see #highlightStrokeOpacity
     * @default {@link JXG.Options.elements.color#strokeColor}
     */
    this.visProp.strokeColor = this.board.options.elements.strokeColor; //'#36393D';

    /**
     * The stroke color of the given geometry element when the user moves the mouse over it.
     * @type string
     * @name JXG.GeometryElement#highlightStrokeColor
     * @see #sstrokeColor
     * @see #strokeWidth
     * @see #strokeOpacity
     * @see #highlightStrokeOpacity
     * @default {@link JXG.Options.elements.color#highlightStrokeColor}
     */
    this.visProp.highlightStrokeColor = this.board.options.elements.highlightStrokeColor;

    /**
     * The fill color of this geometry element.
     * @type string
     * @name JXG.GeometryElement#fillColor
     * @see #highlightFillColor
     * @see #fillOpacity
     * @see #highlightFillOpacity
     * @default {@link JXG.Options.elements.color#fillColor}
     */
    this.visProp.fillColor = this.board.options.elements.fillColor;

    /**
     * The fill color of the given geometry element when the mouse is pointed over it.
     * @type string
     * @name JXG.GeometryElement#highlightFillColor
     * @see #fillColor
     * @see #fillOpacity
     * @see #highlightFillOpacity
     * @default {@link JXG.Options.elements.color#highlightFillColor}
     */
    this.visProp.highlightFillColor = this.board.options.elements.highlightFillColor;

    /**
     * Width of the element's stroke.
     * @type number
     * @name JXG.GeometryElement#strokeWidth
     * @see #strokeColor
     * @see #highlightStrokeColor
     * @see #strokeOpacity
     * @see #highlightStrokeOpacity
     * @default {@link JXG.Options.elements#strokeWidth}
     */
    this.visProp.strokeWidth = this.board.options.elements.strokeWidth;

    /**
     * Opacity for element's stroke color.
     * @type number
     * @name JXG.GeometryElement#strokeOpacity
     * @see #strokeColor
     * @see #highlightStrokeColor
     * @see #strokeWidth
     * @see #highlightStrokeOpacity
     * @default {@link JXG.Options.elements#strokeOpacity}
     */
    this.visProp.strokeOpacity = this.board.options.elements.strokeOpacity;

    /**
     * Opacity for stroke color when the object is highlighted.
     * @type number
     * @name JXG.GeometryElement#highlightStrokeOpacity
     * @see #strokeColor
     * @see #highlightStrokeColor
     * @see #strokeWidth
     * @see #strokeOpacity
     * @default {@link JXG.Options.elements#highlightStrokeOpacity}
     */
    this.visProp.highlightStrokeOpacity = this.board.options.elements.highlightStrokeOpacity;

    /**
     * Opacity for fill color.
     * @type number
     * @name JXG.GeometryElement#fillOpacity
     * @see #fillColor
     * @see #highlightFillColor
     * @see #highlightFillOpacity
     * @default {@link JXG.Options.elements.color#fillOpacity}
     */
    this.visProp.fillOpacity = this.board.options.elements.fillOpacity;

    /**
     * Opacity for fill color when the object is highlighted.
     * @type number
     * @name JXG.GeometryElement#highlightFillOpacity
     * @see #fillColor
     * @see #highlightFillColor
     * @see #fillOpacity
     * @default {@link JXG.Options.elements.color#highlightFillOpacity}
     */
    this.visProp.highlightFillOpacity = this.board.options.elements.highlightFillOpacity;

    /**
     * If true the element will be drawn in grey scale colors to visualize that it's only a draft.
     * @type boolean
     * @name JXG.GeometryElement#draft
     * @default {@link JXG.Options.elements.draft#draft}
     */
    this.visProp.draft = this.board.options.elements.draft.draft;

    /**
     * If false the element won't be visible on the board, otherwise it is shown.
     * @type boolean
     * @name JXG.GeometryElement#visible
     * @see #hideElement
     * @see #showElement
     * @default true
     */
    this.visProp.visible = true;

    /**
     * If true the element will get a shadow.
     * @type boolean
     * @name JXG.GeometryElement#shadow
     * @default false
     */
    this.visProp['shadow'] = false;

    // TODO: withLabel

    // TODO: comment gradient possibilities
    this.visProp['gradient'] = 'none';
    this.visProp['gradientSecondColor'] = 'black';
    this.visProp['gradientAngle'] = '270';
    this.visProp['gradientSecondOpacity'] = this.visProp['fillOpacity'];
    this.visProp['gradientPositionX'] = 0.5;
    this.visProp['gradientPositionY'] = 0.5;
};

/**
 * Add an element as a child to the current element. Can be used to model dependencies between geometry elements.
 * @param {JXG.GeometryElement} obj The dependent object.
 */
JXG.GeometryElement.prototype.addChild = function (obj) {
	var el, el2;

    this.childElements[obj.id] = obj;

    this.addDescendants(obj);

    obj.ancestors[this.id] = this;
    for(el in this.descendants) {
        this.descendants[el].ancestors[this.id] = this;
        for(el2 in this.ancestors) {
            this.descendants[el].ancestors[this.ancestors[el2].id] = this.ancestors[el2];
        }
    }
    for(el in this.ancestors) {
        for(el2 in this.descendants) {
            this.ancestors[el].descendants[this.descendants[el2].id] = this.descendants[el2];
        }
    }
    return this;
};

/**
 * Adds the given object to the descendants list of this object and all its child objects.
 * @param obj The element that is to be added to the descendants list.
 * @private
 * @return
 */
JXG.GeometryElement.prototype.addDescendants = function (/** JXG.GeometryElement */ obj) {
	var el;

    this.descendants[obj.id] = obj;
    for(el in obj.childElements) {
        this.addDescendants(obj.childElements[el]);
    }
    return this;
};

/**
 * Array of strings containing the polynomials defining the element.
 * Used for determining geometric loci the groebner way.
 * @type array
 * @return An array containing polynomials describing the locus of the current object.
 * @private
 */
JXG.GeometryElement.prototype.generatePolynomial = function () {
    return [];
};

/**
 * Animates properties for that object like stroke or fill color, opacity and maybe
 * even more later.
 * @param {Object} hash Object containing propiertes with target values for the animation.
 * @param {number} time Number of milliseconds to complete the animation.
 * @return A reference to the object
 * @type JXG.GeometryElement
 */
JXG.GeometryElement.prototype.animate = function(hash, time) {
    var r, p,
	    delay = 35,
	    steps = Math.ceil(time/(delay * 1.0)),
        i, self = this;

    this.animationData = {};

    var animateColor = function(startRGB, endRGB, property) {
        var hsv1, hsv2, sh, ss, sv;
        hsv1 = JXG.rgb2hsv(startRGB);
        hsv2 = JXG.rgb2hsv(endRGB);

        sh = (hsv2[0]-hsv1[0])/(1.*steps);
        ss = (hsv2[1]-hsv1[1])/(1.*steps);
        sv = (hsv2[2]-hsv1[2])/(1.*steps);
        self.animationData[property] = new Array(steps);
        for(i=0; i<steps; i++) {
            self.animationData[property][steps-i-1] = JXG.hsv2rgb(hsv1[0]+(i+1)*sh, hsv1[1]+(i+1)*ss, hsv1[2]+(i+1)*sv);
        }
    },

    animateFloat = function(start, end, property) {
        start = parseFloat(start);
        end = parseFloat(end);

        // we can't animate without having valid numbers.
        // And parseFloat returns NaN if the given string doesn't contain
        // a valid float number.
        if(isNaN(start) || isNaN(end))
            return;

        var s = (end - start)/(1.*steps);
        self.animationData[property] = new Array(steps);
        for(i=0; i<steps; i++) {
            self.animationData[property][steps-i-1] = start + (i+1)*s;
        }
    };

    for(r in hash) {
        p = r.toLowerCase();
        switch(p) {
            case 'strokecolor':
                    animateColor(this.visProp['strokeColor'], hash[r], 'strokeColor');
                break;
            case 'strokeopacity':
                    animateFloat(this.visProp['strokeOpacity'], hash[r], 'strokeOpacity');
                break;
            case 'strokewidth':
                    animateFloat(this.visProp['strokeWidth'], hash[r], 'strokeWidth');
                break;
            case 'fillcolor':
                    animateColor(this.visProp['fillColor'], hash[r], 'fillColor');
                break;
            case 'fillopacity':
                    animateFloat(this.visProp['fillOpacity'], hash[r], 'fillOpacity');
                break;
        }
    }

	this.board.animationObjects[this.id] = this;
	if(typeof this.board.animationIntervalCode == 'undefined') {
		this.board.animationIntervalCode = window.setInterval('JXG.JSXGraph.boards[\'' + this.board.id + '\'].animate();', delay);
	}

    return this;
};

/**
 * General update method. Should be overwritten by the element itself.
 * Can be used sometimes to commit changes to the object.
 */
JXG.GeometryElement.prototype.update = function() {
    if(this.traced) {
        this.cloneToBackground(true);
    }
    return this;
};

/**
 * Provide updateRenderer method.
 * @private
 */
JXG.GeometryElement.prototype.updateRenderer = function() {
};

/**
 * Hide the element. It will still exist but not visible on the board.
 */
JXG.GeometryElement.prototype.hideElement = function() {
    this.visProp['visible'] = false;
    this.board.renderer.hide(this);
    if (this.label!=null && this.hasLabel) {
        this.label.hiddenByParent = true;
        if(this.label.content.visProp['visible']) {
            this.board.renderer.hide(this.label.content);
        }
    }
    return this;
};

/**
 * Make the element visible.
 */
JXG.GeometryElement.prototype.showElement = function() {
    this.visProp['visible'] = true;
    this.board.renderer.show(this);
    if (this.label!=null && this.hasLabel && this.label.hiddenByParent) {
        this.label.hiddenByParent = false;
        if(this.label.content.visProp['visible']) {
            this.board.renderer.show(this.label.content);
        }
    }
    return this;
};


/* this list is left from the comment below. just to have a list of properties.
* <ul>Possible keys:</ul>
*<li>strokeWidth</li>
*<li>strokeColor</li>
*<li>fillColor</li>
*<li>highlightFillColor</li>
*<li>highlightStrokeColor</li>
*<li>strokeOpacity</li>
*<li>fillOpacity</li>
*<li>highlightFillOpacity</li>
*<li>highlightStrokeOpacity</li>
*<li>labelColor</li>
*<li>visible</li>
*<li>dash</li>
*<li>trace</li>
*<li>style <i>(Point)</i></li>
*<li>fixed</li>
*<li>draft</li>
*<li>showInfobox</li>
*<li>straightFirst <i>(Line)</i></li>
*<li>straightLast <i>(Line)</i></li>
*<li>firstArrow <i>(Line,Arc)</li>
*<li>lastArrow <i>(Line,Arc)</li>
*<li>withTicks <i>(Line)</li>
*</ul>*/

/**
 * Sets an arbitrary number of properties.
 * @param % Arbitrary number of strings, containing "key:value" pairs.
 * The possible key values are the element and class fields in this documentation.
 * @example
 * // Set property directly on creation of an element using the attributes object parameter
 * var board = JXG.JSXGraph.initBoard('jxgbox', {boundingbox: [-1, 5, 5, 1]};
 * var p = board.createElement('point', [2, 2], {visible: false});
 *
 * // Now make this point visible and fixed:
 * p.setProperty('fixed:true', 'visible:true');
 *
 * // Alternatively you can use #hideElement resp. #showElement:
 * p.hideElement();
 */
JXG.GeometryElement.prototype.setProperty = function () {
    var i, key, color, pairRaw,
        opacity,
        pair;

    for (i=0; i<arguments.length; i++) {
        pairRaw = arguments[i];
        if (typeof pairRaw == 'string') {    // pairRaw is string of the form 'key:value'
            pair = pairRaw.split(':');
            // trim pair[0] and pair[1]
            pair[0] = pair[0].replace (/^\s+/, '').replace (/\s+$/, '');
            pair[1] = pair[1].replace (/^\s+/, '').replace (/\s+$/, '');
        } else if (!JXG.isArray(pairRaw)) {    // pairRaw consists of objects of the form {key1:value1,key2:value2,...}
            /*
            for (var i=0; i<Object.keys(pairRaw).length;i++) {  // Here, the prototype lib is used (Object.keys, Object.isArray)
                var key = Object.keys(pairRaw)[i];
                this.setProperty([key,pairRaw[key]]);
            }
            */
            for (key in pairRaw) {
                this.setProperty([key,pairRaw[key]]);
            }
            return this;
        } else {                             // pairRaw consists of array [key,value]
            pair = pairRaw;
        }
        if (pair[1]==null) continue;
        switch(pair[0].replace(/\s+/g).toLowerCase()) {   // Whitespace entfernt und in Kleinbuchstaben umgewandelt.
            case 'strokewidth':
                this.visProp['strokeWidth'] = pair[1];
                this.board.renderer.setObjectStrokeWidth(this, this.visProp['strokeWidth']);
                break;
            case 'strokecolor':
                color = pair[1];
                if (color.length=='9' && color.substr(0,1)=='#') {
                    opacity = color.substr(7,2);
                    color = color.substr(0,7);
                }
                else {
                    opacity = 'FF';
                }
                this.visProp['strokeColor'] = color;
                this.visProp['strokeOpacity'] = parseInt(opacity.toUpperCase(),16)/255;
                this.board.renderer.setObjectStrokeColor(this, this.visProp['strokeColor'], this.visProp['strokeOpacity']);
                break;
            case 'fillcolor':
                color = pair[1];
                if (color.length=='9' && color.substr(0,1)=='#') {
                    opacity = color.substr(7,2);
                    color = color.substr(0,7);
                }
                else {
                    opacity = 'FF';
                }
                this.visProp['fillColor'] = color;
                this.visProp['fillOpacity'] = parseInt(opacity.toUpperCase(),16)/255;
                this.board.renderer.setObjectFillColor(this, this.visProp['fillColor'], this.visProp['fillOpacity']);
                break;
            case 'highlightstrokecolor':
                color = pair[1];
                if (color.length=='9' && color.substr(0,1)=='#') {
                    opacity = color.substr(7,2);
                    color = color.substr(0,7);
                }
                else {
                    opacity = 'FF';
                }
                this.visProp['highlightStrokeColor'] = color;
                this.visProp['highlightStrokeOpacity'] = parseInt(opacity.toUpperCase(),16)/255;
                break;
            case 'highlightfillcolor':
                color = pair[1];
                if (color.length=='9' && color.substr(0,1)=='#') {
                    opacity = color.substr(7,2);
                    color = color.substr(0,7);
                }
                else {
                    opacity = 'FF';
                }
                this.visProp['highlightFillColor'] = color;
                this.visProp['highlightFillOpacity'] = parseInt(opacity.toUpperCase(),16)/255;
                break;
            case 'fillopacity':
                this.visProp['fillOpacity'] = pair[1];
                this.board.renderer.setObjectFillColor(this, this.visProp['fillColor'], this.visProp['fillOpacity']);
                break;
            case 'strokeopacity':
                this.visProp['strokeOpacity'] = pair[1];
                this.board.renderer.setObjectStrokeColor(this, this.visProp['strokeColor'], this.visProp['strokeOpacity']);
                break;
            case 'highlightfillopacity':
                this.visProp['highlightFillOpacity'] = pair[1];
                break;
            case 'highlightstrokeopacity':
                this.visProp['highlightStrokeOpacity'] = pair[1];
                break;
            case 'labelcolor':
                color = pair[1];
                if (color.length=='9' && color.substr(0,1)=='#') {
                    opacity = color.substr(7,2);
                    color = color.substr(0,7);
                }
                else {
                    opacity = 'FF';
                }
                if(opacity == '00') {
                    if (this.label!=null && this.hasLabel) {
                        this.label.content.hideElement();
                    }
                }
                if(this.label!=null && this.hasLabel) {
                    this.label.color = color;
                    this.board.renderer.setObjectStrokeColor(this.label.content, color, opacity);
                }
                if(this.type == JXG.OBJECT_TYPE_TEXT) {
                    this.visProp['strokeColor'] = color;
                    this.board.renderer.setObjectStrokeColor(this, this.visProp['strokeColor'], 1);
                }
                break;
            case 'showinfobox':
                if(pair[1] == 'false' || pair[1] == false) {
                    this.showInfobox = false;
                }
                else if(pair[1] == 'true' || pair[1] == true) {
                    this.showInfobox = true;
                }
                break;
            case 'visible':
                if(pair[1] == 'false' || pair[1] == false) {
                    this.visProp['visible'] = false;
                    this.hideElement();
                }
                else if(pair[1] == 'true' || pair[1] == true) {
                    this.visProp['visible'] = true;
                    this.showElement();
                }
                break;
            case 'dash':
                this.setDash(pair[1]);
                break;
            case 'trace':
                if(pair[1] == 'false' || pair[1] == false) {
                    this.traced = false;
                }
                else if(pair[1] == 'true' || pair[1] == true) {
                    this.traced = true;
                }
                break;
            case 'style':
                this.setStyle(1*pair[1]);
                break;
            case 'face':
            	if(this.elementClass == JXG.OBJECT_CLASS_POINT)
            		this.setFace(pair[1]);
                break;
            case 'size':
            	if(this.elementClass == JXG.OBJECT_CLASS_POINT) {
            		this.visProp['size'] = 1*pair[1];
                	this.board.renderer.updatePoint(this);
        		}
                break;
            case 'fixed':
                this.fixed = ((pair[1]=='false') || (pair[1]==false)) ? false : true;
                break;
            case 'shadow':
                if(pair[1] == 'false' || pair[1] == false) {
                    this.visProp['shadow'] = false;
                }
                else if(pair[1] == 'true' || pair[1] == true) {
                    this.visProp['shadow'] = true;
                }
                this.board.renderer.setShadow(this);
                break;
            case 'gradient':
                this.visProp['gradient'] = pair[1];
                this.board.renderer.setGradient(this);
                break;
            case 'gradientsecondcolor':
                color = pair[1];
                if (color.length=='9' && color.substr(0,1)=='#') {
                    opacity = color.substr(7,2);
                    color = color.substr(0,7);
                }
                else {
                    opacity = 'FF';
                }
                this.visProp['gradientSecondColor'] = color;
                this.visProp['gradientSecondOpacity'] = parseInt(opacity.toUpperCase(),16)/255;
                this.board.renderer.updateGradient(this);
                break;
            case 'gradientsecondopacity':
                this.visProp['gradientSecondOpacity'] = pair[1];
                this.board.renderer.updateGradient(this);
                break;
            case 'draft':
                if(pair[1] == 'false' || pair[1] == false) {
                    if(this.visProp['draft'] == true) {
                        this.visProp['draft'] = false;
                        this.board.renderer.removeDraft(this);
                    }
                }
                else if(pair[1] == 'true' || pair[1] == true) {
                    this.visProp['draft'] = true;
                    this.board.renderer.setDraft(this);
                }
                break;
            case 'straightfirst':
                if(pair[1] == 'false' || pair[1] == false) {
                    this.visProp['straightFirst'] = false;
                }
                else if(pair[1] == 'true' || pair[1] == true) {
                    this.visProp['straightFirst'] = true;
                }
                this.setStraight(this.visProp['straightFirst'], this.visProp['straightLast']);
                break;
            case 'straightlast':
                if(pair[1] == 'false' || pair[1] == false) {
                    this.visProp['straightLast'] = false;
                }
                else if(pair[1] == 'true' || pair[1] == true) {
                    this.visProp['straightLast'] = true;
                }
                this.setStraight(this.visProp['straightFirst'], this.visProp['straightLast']);
                break;
            case 'firstarrow':
                if(pair[1] == 'false' || pair[1] == false) {
                    this.visProp['firstArrow'] = false;
                }
                else if(pair[1] == 'true' || pair[1] == true) {
                    this.visProp['firstArrow'] = true;
                }
                this.setArrow(this.visProp['firstArrow'], this.visProp['lastArrow']);
                break;
            case 'lastarrow':
                if(pair[1] == 'false' || pair[1] == false) {
                    this.visProp['lastArrow'] = false;
                }
                else if(pair[1] == 'true' || pair[1] == true) {
                    this.visProp['lastArrow'] = true;
                }
                this.setArrow(this.visProp['firstArrow'], this.visProp['lastArrow']);
                break;
            case 'curvetype':
                this.curveType = pair[1];
                break;
            case 'fontsize':
                this.visProp['fontSize'] = pair[1];
                break;
            case 'insertticks':
                if(this.type == JXG.OBJECT_TYPE_TICKS) {
                    var old = this.insertTicks;
                    this.insertTicks = true;
                    if(pair[1] == 'false' || pair[1] == false) {
                        this.insertTicks = false;
                    }
                    if(old != this.insertTicks) this.calculateTicksCoordinates();
                }
                break;
            case 'drawlabels':
                if(this.type == JXG.OBJECT_TYPE_TICKS) {
                    var old = this.drawLabels;
                    this.drawLabels = true;
                    if(pair[1] == 'false' || pair[1] == false) {
                        this.drawLabels = false;
                    }
                    if(old != this.drawLabels) this.calculateTicksCoordinates();
                }
                break;
            case 'drawzero':
                if(this.type == JXG.OBJECT_TYPE_TICKS) {
                    var old = this.drawZero;
                    this.drawZero = true;
                    if(pair[1] == 'false' || pair[1] == false) {
                        this.drawZero = false;
                    }
                    if(old != this.drawZero) this.calculateTicksCoordinates();
                }
                break;
            case 'minorticks':
                if(this.type == JXG.OBJECT_TYPE_TICKS) {
                    var old = this.minorTicks;
                    if((pair[1] != null) && (pair[1] > 0))
                        this.minorTicks = pair[1];
                    if(old != this.minorTicks) this.calculateTicksCoordinates();
                }
                break;
            case 'majortickheight':
                if(this.type == JXG.OBJECT_TYPE_TICKS) {
                    var old = this.majorHeight;
                    if((pair[1] != null) && (pair[1] > 0))
                        this.majorHeight = pair[1];
                    if(old != this.majorHeight) this.calculateTicksCoordinates();
                }
                break;
            case 'minortickheight':
                if(this.type == JXG.OBJECT_TYPE_TICKS) {
                    var old = this.minorHeight;
                    if((pair[1] != null) && (pair[1] > 0))
                        this.minorHeight = pair[1];
                    if(old != this.minorHeight) this.calculateTicksCoordinates();
                }
                break;
            case 'snapwidth':
                if(this.type == JXG.OBJECT_TYPE_GLIDER) {
                    this.snapWidth = pair[1];
                }
        }
    }
    return this;
};

/**
 * Set the dash style of an object. See {@link #dash} for a list of available dash styles.
 * You should use {@link #setProperty} instead of this method.
 * @param {number} dash Indicates the new dash style
 * @private
*/
JXG.GeometryElement.prototype.setDash = function(dash) {
    this.visProp['dash'] = dash;
    this.board.renderer.setDashStyle(this,this.visProp);
    return this;
};

/**
 * Notify all child elements for updates.
 * @private
 */
JXG.GeometryElement.prototype.prepareUpdate = function() {
    this.needsUpdate = true;
    return this; // Im Moment steigen wir nicht rekursiv hinab
    /* End of function  */

    /*
    var el;
    for(el in this.childElements) {
        // Wurde das Element vielleicht geloescht?
        if(this.board.objects[el] != undefined) {
            // Nein, wurde es nicht, also updaten
            this.childElements[el].prepareUpdate();
        } else { //  es wurde geloescht, also aus dem Array entfernen
            delete(this.childElements[el]);
        }
    }
    */
};

/**
 * Removes the element from the construction.
 */
JXG.GeometryElement.prototype.remove = function() {
    this.board.renderer.remove(document.getElementById(this.id));
    if (this.hasLabel) {
        this.board.renderer.remove(document.getElementById(this.label.content.id));
    }
    return this;
};

/**
 * Returns the coords object where a text that is bound to the element shall be drawn.
 * Differs in some cases from the values that getLabelAnchor returns.
 * @type JXG.Coords
 * @return JXG.Coords Place where the text shall be drawn.
 * @see #getLabelAnchor
 * @private
 */
JXG.GeometryElement.prototype.getTextAnchor = function() {
    return new JXG.Coords(JXG.COORDS_BY_USER, [0,0], this.board);
};

/**
 * Returns the coords object where the label of the element shall be drawn.
  * Differs in some cases from the values that getTextAnchor returns.
 * @type JXG.Coords
 * @return JXG.Coords Place where the label of an element shall be drawn.
  * @see #getTextAnchor
 * @private
 */
JXG.GeometryElement.prototype.getLabelAnchor = function() {
    return new JXG.Coords(JXG.COORDS_BY_USER, [0,0], this.board);
};

/**
 * TODO
 * Was hat das hier verloren? Styles gibts doch nur fuer Punkte oder?
 * Sollte das dann nicht nur in Point.js zu finden sein? --michael
 * @private
 */
JXG.GeometryElement.prototype.setStyle = function(x) {
    return this;
};

/**
 * TODO
 * Was hat das hier verloren? "Straights" gibts doch nur fuer Lines oder?
 * Sollte das dann nicht nur in Line.js zu finden sein? --michael
 * @private
 */
JXG.GeometryElement.prototype.setStraight = function(x,y) {
    return this;
};

/**
 * TODO
 * Dito setStraight. Das gilt doch eh nur fuer lines, also wozu hier reinstellen? --michael
 * @private
 */
JXG.GeometryElement.prototype.setArrow = function(firstArrow,lastArrow) {
    return this;
};

/**
 * Creates a label element for this geometry element.
 * Doesn't add the label to the board, so it shouldn't be called itself. Use {@link #addLabelToElement} instead.
 * @param {boolean} withLabel true if a label shall be initialized, false otherwise.
 * @see #addLabelToElement
 * @private
 */
JXG.GeometryElement.prototype.createLabel = function(withLabel,coords) {
    var isTmpId = false;
    if (typeof coords=='undefined' || coords==null) {
        coords = [10,10];
    }
    this.nameHTML = this.board.algebra.replaceSup(this.board.algebra.replaceSub(this.name));
    this.label = {};
    if (typeof withLabel=='undefined' || withLabel==true) {
        if (this.board.objects[this.id]==null) {
            this.board.objects[this.id] = this;
            isTmpId = true;
        }
        this.label.relativeCoords = coords;
        this.label.content = new JXG.Text(this.board, this.nameHTML, this.id,
            [this.label.relativeCoords[0]/(this.board.stretchX),this.label.relativeCoords[1]/(this.board.stretchY)], this.id+"Label", "", null, true, this.board.options.text.defaultType);
        if (isTmpId) delete(this.board.objects[this.id]);
        this.label.color = '#000000';
        if(!this.visProp['visible']) {
            this.label.hiddenByParent = true;
            this.label.content.visProp['visible'] = false;
        }
        this.hasLabel = true;
    }
    return this;
};

/**
 * Adds a label to the element.
 */
JXG.GeometryElement.prototype.addLabelToElement = function() {
    this.createLabel(true);
    this.label.content.id = this.id+"Label";
    this.board.addText(this.label.content);
    this.board.renderer.drawText(this.label.content);
    if(!this.label.content.visProp['visible']) {
        board.renderer.hide(this.label.content);
    }
    return this;
};

/**
 * Highlights the element.
 */
JXG.GeometryElement.prototype.highlight = function() {
    this.board.renderer.highlight(this);
    return this;
};

/**
 * Uses the "normal" properties of the element.
 */
JXG.GeometryElement.prototype.noHighlight = function() {
    this.board.renderer.noHighlight(this);
    return this;
};

/**
 * Removes all objects generated by the trace function.
 */
JXG.GeometryElement.prototype.clearTrace = function() {
    var obj;

    for(obj in this.traces) {
        this.board.renderer.remove(this.traces[obj]);
    }
    this.numTraces = 0;
    return this;
};

/**
 * Copy element to background. Has to be implemented in the element itself.
 * @private
 */
JXG.GeometryElement.prototype.cloneToBackground = function(addToTrace) {
    return this;
};

// [c,b0,b1,a,k]
/**
 * Normalize the element's standard form.
 * @private
 */
JXG.GeometryElement.prototype.normalize = function() {
    this.stdform = this.board.algebra.normalize(this.stdform);
    return this;
};

/**
 * EXPERIMENTAL. Generate JSON object code of visProp and other properties.
 * @type string
 * @private
 * @return JSON string containing element's properties.
 */
JXG.GeometryElement.prototype.toJSON = function() {
    var json = '{"name":' + this.name;
    json += ', ' + '"id":' + this.id;

    var vis = [];
    for (var key in this.visProp) {
        if (this.visProp[key]!=null) {
            vis.push('"' + key + '":' + this.visProp[key]);
        }
    }
    json += ', "visProp":{'+vis.toString()+'}';
    json +='}';

    return json;
};

/**
  * Setting visPropOld is done in an none object oriented version
  * since otherwise there would be problems in cloneToBackground
  */
JXG.clearVisPropOld = function(el) {
    el.visPropOld = {};
    el.visPropOld['strokeColor']= '';
    el.visPropOld['strokeOpacity']= '';
    el.visPropOld['strokeWidth']= '';
    el.visPropOld['fillColor']= '';
    el.visPropOld['fillOpacity']= '';
    el.visPropOld['shadow']= false;
    el.visPropOld['firstArrow'] = false;
    el.visPropOld['lastArrow'] = false;
};

/* 
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/** 
 * @fileoverview In this file the Coords object is defined, a class to manage all
 * properties and methods coordinates usually have.
 * @author graphjs
 * @version 0.1
 */

JXG.COORDS_BY_USER = 0x0001;
JXG.COORDS_BY_SCREEN = 0x0002;

/**
 * Constructs a new Coordinates object.
 * @class This is the Coordinates class.  
 * All members a coordinate has to provide
 * are defined here.
 * @param {int} method The type of coordinates given by the user. Accepted values are <b>COORDS_BY_SCREEN</b> and <b>COORDS_BY_USER</b>.
 * @param {Array} coordinates An array of affine coordinates.
 * @param {JXG.AbstractRenderer} renderer A reference to a Renderer.
 * @constructor
 */
JXG.Coords = function (method, coordinates, board) {
    /**
     * Stores the board the object is used on.
     * @type JXG.Board
     */
    this.board = board;
    
    /**
     * Stores coordinates for user view as homogeneous coordinates.
     * @type Array
     */
    this.usrCoords = [];
    /**
     * Stores coordinates for screen view as homogeneous coordinates.
     * @type Array
     */
    this.scrCoords = [];
    
    if(method == JXG.COORDS_BY_USER) {
        if (coordinates.length<=2) {
            this.usrCoords[0] = 1.0;
            this.usrCoords[1] = coordinates[0];
            this.usrCoords[2] = coordinates[1];
        } else {  // homogeneous coordinates
            this.usrCoords[0] = coordinates[0];
            this.usrCoords[1] = coordinates[1];
            this.usrCoords[2] = coordinates[2];
            this.normalizeUsrCoords();
        }
        this.usr2screen();
    } else {
        this.scrCoords[0] = 1.0;
        this.scrCoords[1] = coordinates[0];
        this.scrCoords[2] = coordinates[1];
        this.screen2usr();
    }
};

/**
    * Normalize homogeneous coordinates
    * @private
    */
JXG.Coords.prototype.normalizeUsrCoords = function() {
    var eps = 0.000001;
    if (Math.abs(this.usrCoords[0])>eps) {
        this.usrCoords[1] /= this.usrCoords[0];
        this.usrCoords[2] /= this.usrCoords[0];
        this.usrCoords[0] = 1.0;
    }
};

/**
 * Compute screen coordinates out of given user coordinates.
 * @private
 */
JXG.Coords.prototype.usr2screen = function(doRound) {
    var mround = Math.round,  // Is faster on IE, maybe slower with JIT compilers
        b = this.board,
        uc = this.usrCoords,
        oc = this.board.origin.scrCoords;
        
    if (doRound==null || doRound) {
        this.scrCoords[0] = mround(uc[0]);
        this.scrCoords[1] = mround(uc[0]*oc[1] + uc[1]*b.stretchX);
        this.scrCoords[2] = mround(uc[0]*oc[2] - uc[2]*b.stretchY);
    } else {
        this.scrCoords[0] = uc[0];
        this.scrCoords[1] = uc[0]*oc[1] + uc[1]*b.stretchX;
        this.scrCoords[2] = uc[0]*oc[2] - uc[2]*b.stretchY;
    }
};

/**
 * Compute user coordinates out of given screen coordinates.
 * @private
 */
JXG.Coords.prototype.screen2usr = function() {
    var o = this.board.origin.scrCoords,
        sc = this.scrCoords,
        b = this.board;
    this.usrCoords[0] =  1.0;
    this.usrCoords[1] = (sc[1] - o[1])/b.stretchX;
    this.usrCoords[2] = (o[2] - sc[2])/b.stretchY;
};

/**
 * Calculate distance of one point to another.
 * @param {int} method The type of coordinates used here. Possible values are <b>JXG.COORDS_BY_USER</b> and <b>JXG.COORDS_BY_SCREEN</b>.
 * @param {JXG.Coords} coordinates The Coords object to which the distance is calculated.
 */
JXG.Coords.prototype.distance = function(meth, crd) {
    var sum = 0,
        c,
        ucr = this.usrCoords,
        scr = this.scrCoords,
        f;
        
    if (meth == JXG.COORDS_BY_USER) {
        c = crd.usrCoords;
        f = ucr[0]-c[0];
        sum = f*f;
        f = ucr[1]-c[1];
        sum += f*f;
        f = ucr[2]-c[2];
        sum += f*f;
    } else {
        c = crd.scrCoords;
        f = scr[0]-c[0];
        sum = f*f;
        f = scr[1]-c[1];
        sum += f*f;
        f = scr[2]-c[2];
        sum += f*f;
    }

    /*
    if (meth == JXG.COORDS_BY_USER) {
//        if (Math.abs(this.usrCoords[0]+coordinates.usrCoords[0])>eps) {
//            return Infinity;
//        }
        for (i=0; i<=this.board.dimension; i++) {
            f = this.usrCoords[i] - coordinates.usrCoords[i];
            sum += f*f;
        }
    } else {
//        if (Math.abs(this.scrCoords[0]+coordinates.scrCoords[0])>eps) {
//            return Infinity;
//        }
        for (i=0; i<=this.board.dimension; i++) {
            f = this.scrCoords[i] - coordinates.scrCoords[i];
            sum += f*f;
        }
    }
    */
    
    return Math.sqrt(sum);
};

/**
 * Set coordinates by method
 * @param {int} method The type of coordinates used here. Possible values are <b>COORDS_BY_USER</b> and <b>COORDS_BY_SCREEN</b>.
 * @param {Array} coordinates An array of affine coordinates the Coords object is set to.
 * @param {boolean} optional flag If true or null round the coordinates in usr2screen. This is used in smooth curve plotting.
 * The IE needs rounded coordinates. Id doRound==false we have to round in updatePathString.
 */
JXG.Coords.prototype.setCoordinates = function(method, crd, doRound) {
    var uc = this.usrCoords,
        sc = this.scrCoords;
        
    if (method == JXG.COORDS_BY_USER) {
        if (crd.length==2) { // Euclidean coordinates
            uc[0] = 1.0;
            uc[1] = crd[0];
            uc[2] = crd[1];
        } else { // Homogeneous coordinates (normalized)
            uc[0] = crd[0];
            uc[1] = crd[1];
            uc[2] = crd[2];
            this.normalizeUsrCoords();
        }
        this.usr2screen(doRound);
    } else {
        sc[1] = crd[0];
        sc[2] = crd[1];
        this.screen2usr();
    }
};

/*
    Copyright 2008-2010
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview The geometry object Point is defined in this file. Point stores all
 * style and functional properties that are required to draw and move a point on
 * a board.
 * @author graphjs
 * @version 0.1
 */


JXG.POINT_STYLE_X_SMALL      = 0;  // a small sized x
JXG.POINT_STYLE_X            = 1;  // a medium sized x
JXG.POINT_STYLE_X_BIG        = 2;  // a big sized x
JXG.POINT_STYLE_CIRCLE_TINY  = 3;  // a tiny circle
JXG.POINT_STYLE_CIRCLE_SMALL = 4;  // a small circle
JXG.POINT_STYLE_CIRCLE       = 5;  // a medium circle
JXG.POINT_STYLE_CIRCLE_BIG   = 6;  // a big circle
JXG.POINT_STYLE_SQUARE_SMALL = 7;  // a small rectangle
JXG.POINT_STYLE_SQUARE       = 8;  // a medium rectangle
JXG.POINT_STYLE_SQUARE_BIG   = 9;  // a big rectangle
JXG.POINT_STYLE_PLUS_SMALL   = 10; // a small +
JXG.POINT_STYLE_PLUS         = 11; // a medium +
JXG.POINT_STYLE_PLUS_BIG     = 12; // a big +

/**
 * A point is the basic geometric element. Based on points lines and circles can be constructed which can be intersected
 * which in turn are points again which can be used to construct new lines, circles, polygons, etc. This class holds methods for
 * all kind of points like free points, gliders, and intersection points.
 * @class Creates a new point object. Do not use this constructor to create a point. Use {@link JXG.Board#create} with
 * type {@link Point}, {@link Glider}, or {@link Intersection} instead.  
 * @augments JXG.GeometryElement
 * @param {string,JXG.Board} board The board the new point is drawn on.
 * @param {Array} coordinates An array with the affine user coordinates of the point.
 * @param {String} id Unique identifier for the point. If null or an empty string is given,
 *  an unique id will be generated by Board
 * @param {String} name Not necessarily unique name for the point. If null or an
 *  empty string is given, an unique name will be generated
 * @param {boolean} show False if the point is invisible, True otherwise
 * @see JXG.Board#generateName
 * @see JXG.Board#addPoint
 */
JXG.Point = function (board, coordinates, id, name, show, withLabel, layer) {
    this.constructor();
    
    /**
     * Type of point; Possible values are {@link JXG.OBJECT_TYPE_POINT}, {@link JXG.OBJECT_TYPE_GLIDER}, {@link JXG.OBJECT_TYPE_CAS}.
     * @default {@link JXG.OBJECT_TYPE_POINT}
     * @type number
     * @private
     */
    this.type = JXG.OBJECT_TYPE_POINT;
    
    /**
     * Class of this point element; Values is OBJECT_CLASS_POINT.
     * @constant
     * @type number
     * @private
     */
    this.elementClass = JXG.OBJECT_CLASS_POINT;

    this.init(board, id, name);

    if (coordinates==null) {
        coordinates=[0,0];
    }
    /**
     * Coordinates of the point.
     * @type JXG.Coords
     * @private
     */
    this.coords = new JXG.Coords(JXG.COORDS_BY_USER, coordinates, this.board);
    this.initialCoords = new JXG.Coords(JXG.COORDS_BY_USER, coordinates, this.board);

    /**
     * Set the display layer.
     */
    if (layer == null) layer = board.options.layer['point'];
    this.layer = layer;

    /**
     * If true, the infobox is shown on mouse over, else not.
     * @type boolean
     * @default true
     */
    this.showInfobox = true;
    
    /**
     * Descriptive character, displayed next to the point
     * @type JXG.Label
     * @private
     */
    this.label = {};
    this.label.relativeCoords = [10,10];
    this.nameHTML = this.board.algebra.replaceSup(this.board.algebra.replaceSub(this.name)); //?
    if (typeof withLabel=='undefined' || withLabel==true) {
        this.board.objects[this.id] = this;
        this.label.content = new JXG.Text(this.board, this.nameHTML, this.id, 
            [this.label.relativeCoords[0]/this.board.stretchX,this.label.relativeCoords[1]/this.board.stretchY], this.id+"Label", "", null, true, this.board.options.text.defaultType);
        delete(this.board.objects[this.id]);

        this.label.color = '#000000';
        if(!show) {
            this.label.hiddenByParent = true;
            this.label.content.visProp['visible'] = false;
        }
        this.hasLabel = true;
    } else {
        this.showInfobox = false;
    }
    
    /**
     * False: Point can be moved, True: Point can't be moved with the mouse.
     * @type boolean
     * @default false
     */
    this.fixed = false;
    
    /**
     * Relative position on a line if point is a glider on a line.
     * @type number
     * @private
     */
    this.position = null;

    /**
     * Determines whether the point slides on a polygon if point is a glider.
     * @type boolean
     * @default false
     * @private
     */
    this.onPolygon = false;
    
    /**
     * There are different point styles which differ in appearance and size.
     * Possible values are
     * <table><tr><th>Constant name</th><th>Value</th><th>Meaning</th><th>face</th><th>size</th></tr>
     * <tr><td>JXG.POINT_STYLE_X_SMALL</td><td>0</td><td>small sized x</td><td>cross or x</td><td>2</td></tr>
     * <tr><td>JXG.POINT_STYLE_X</td><td>1</td><td>medium sized x</td><td>cross or x</td><td>3</td></tr>
     * <tr><td>JXG.POINT_STYLE_X_BIG</td><td>2</td><td>big sized x</td><td>cross or x</td><td>4</td></tr>
     * <tr><td>JXG.POINT_STYLE_CIRCLE_TINY</td><td>3</td><td>tiny circle</td><td>circle or o</td><td>1</td></tr>
     * <tr><td>JXG.POINT_STYLE_CIRCLE_SMALL</td><td>4</td><td>small circle</td><td>circle or o</td><td>2</td></tr>
     * <tr><td>JXG.POINT_STYLE_CIRCLE</td><td>5</td><td>medium circle</td><td>circle or o</td><td>3</td></tr>
     * <tr><td>JXG.POINT_STYLE_CIRCLE_BIG</td><td>6</td><td>big circle</td><td>circle or o</td><td>4</td></tr>
     * <tr><td>JXG.POINT_STYLE_SQUARE_SMALL</td><td>7</td><td>small rectangle</td><td>square or []</td><td>2</td></tr>
     * <tr><td>JXG.POINT_STYLE_SQUARE</td><td>8</td><td>medium rectangle</td><td>square or []</td><td>3</td></tr>
     * <tr><td>JXG.POINT_STYLE_SQUARE_BIG</td><td>9</td><td>big rectangle</td><td>square or []</td><td>4</td></tr>
     * <tr><td>JXG.POINT_STYLE_PLUS_SMALL</td><td>10</td><td>small +</td><td>plus or +</td><td>2</td></tr>
     * <tr><td>JXG.POINT_STYLE_PLUS</td><td>11</td><td>medium +</td><td>plus or +</td><td>3</td></tr>
     * <tr><td>JXG.POINT_STYLE_PLUS_BIG</td><td>12</td><td>big +</td><td>plus or +</td><td>4</td></tr></table>
     * <b>Hint:</b> This attribute is internally replaced by face and size, whose opportunities are wider, , as given in the table above.
     * @see JXG.Point#face
     * @see JXG.Point#size
     * @type number
     * @see #setStyle
     * @default JXG.Options.point#style
     * @name JXG.Point#style
     * @deprecated
     */
    this.visProp['style'] = this.board.options.point.style;
    
    /**
     * There are different point styles which differ in appearance.
     * Posssible values are
     * <table><tr><th>Value</th></tr>
     * <tr><td>cross</td></tr>
     * <tr><td>circle</td></tr>
     * <tr><td>square</td></tr>
     * <tr><td>plus</td></tr>
     * <tr><td>diamond</td></tr>
     * <tr><td>triangleUp</td></tr>
     * <tr><td>triangleDown</td></tr>
     * <tr><td>triangleLeft</td></tr>
     * <tr><td>triangleRight</td></tr>
     * </table>
     * @type string
     * @see #setStyle
     * @default circle
     * @name JXG.Point#face
     */
    this.visProp['face'] = 'circle';
    /**
    * Determines the size of a point.
    * Means radius resp. half the width of a point (depending on the face).
     * @see JXG.Point#face    
    * @type number
     * @see #setStyle
     * @default 3     
     * @name JXG.Point#size
     */    
    this.visProp['size'] = 3;

    /**
     * Size of the point. This is just for the renderer and the hasPoint() method
     * to draw the point as a circle.
     * @type number
     * @private
     */
    //this.r = this.board.options.precision.hasPoint;
    
    /*
     * The visprop properties are documented in JXG.GeometryElement
     */
    this.visProp['fillColor'] = this.board.options.point.fillColor;
    this.visProp['highlightFillColor'] = this.board.options.point.highlightFillColor;  
    this.visProp['strokeColor'] = this.board.options.point.strokeColor;
    this.visProp['highlightStrokeColor'] = this.board.options.point.highlightStrokeColor; 
    this.visProp['strokeWidth'] = this.board.options.point.strokeWidth;    
    this.visProp['visible'] = show; 

    /**
     * When used as a glider this member stores the object, where to glide on. To set the object to glide on use the method
     * {@link JXG.Point#makeGlider} and DO NOT set this property directly as it will break the dependency tree.
     * TODO: Requires renaming to glideObject
     * @type JXG.GeometryElement
     * @name Glider#slideObject
     */
    this.slideObject = null;
    
    /**
     * Stores the groups of this point in an array of Group.
     * @type array
     * @see JXG.Group
     * @private
     */
    this.group = [];
    
    /* Register point at board. */
    this.id = this.board.addPoint(this);
};

/**
 * Inherits here from {@link JXG.GeometryElement}.
 */
JXG.Point.prototype = new JXG.GeometryElement();

/**
 * Checks whether (x,y) is near the point.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @type boolean
 * @return True if (x,y) is near the point, False otherwise.
 * @private
 */
JXG.Point.prototype.hasPoint = function (x,y) {
    var coordsScr = this.coords.scrCoords, r;
    r = this.visProp['size'];
    if(r < this.board.options.precision.hasPoint) {
        r = this.board.options.precision.hasPoint;
    }
    return ((Math.abs(coordsScr[1]-x) < r+2) && (Math.abs(coordsScr[2]-y)) < r+2);
};

/**
* Dummy function for unconstrained points or gliders.
* @private
*/
JXG.Point.prototype.updateConstraint = function() { return this; };

/**
 * Updates the position of the point.
 */
JXG.Point.prototype.update = function (fromParent) {
    if (!this.needsUpdate) { return; }

    if(typeof fromParent == 'undefined') {
        fromParent = false;
    }
  
    if(this.traced) {
        this.cloneToBackground(true);
    }
 /*
     * We need to calculate the new coordinates no matter of the points visibility because
     * a child could be visible and depend on the coordinates of the point (e.g. perpendicular).
     * 
     * Check if point is a glider and calculate new coords in dependency of this.slideObject.
     * This function is called with fromParent==true for example if
     * the defining elements of the line or circle have been changed.
     */
    if(this.type == JXG.OBJECT_TYPE_GLIDER) {
        if(this.slideObject.type == JXG.OBJECT_TYPE_CIRCLE) {
//fromParent = false;        
            if (fromParent) {
                this.coords.setCoordinates(JXG.COORDS_BY_USER, [this.slideObject.midpoint.X()+Math.cos(this.position),this.slideObject.midpoint.Y()+Math.sin(this.position)]);
                this.coords  = this.board.algebra.projectPointToCircle(this, this.slideObject);
            } else {
                this.coords  = this.board.algebra.projectPointToCircle(this, this.slideObject);
                this.position = this.board.algebra.rad([this.slideObject.midpoint.X()+1.0,this.slideObject.midpoint.Y()],this.slideObject.midpoint,this);
            }
        } else if(this.slideObject.type == JXG.OBJECT_TYPE_LINE) {
            this.coords  = this.board.algebra.projectPointToLine(this, this.slideObject);
            
            var p1coords = this.slideObject.point1.coords;
            var p2coords = this.slideObject.point2.coords;
            if (fromParent) {
                if (Math.abs(p1coords.usrCoords[0])>=JXG.Math.eps && Math.abs(p2coords.usrCoords[0])>=JXG.Math.eps) {
                    this.coords.setCoordinates(JXG.COORDS_BY_USER, 
                                           [p1coords.usrCoords[1] + this.position*(p2coords.usrCoords[1] - p1coords.usrCoords[1]),
                                            p1coords.usrCoords[2] + this.position*(p2coords.usrCoords[2] - p1coords.usrCoords[2])]);
                }
            } else {
                var factor = 1;
                var distP1S = p1coords.distance(JXG.COORDS_BY_USER, this.coords);
                var distP1P2 = p1coords.distance(JXG.COORDS_BY_USER, p2coords);
                var distP2S = p2coords.distance(JXG.COORDS_BY_USER, this.coords);
                
                if( ((distP1S > distP1P2) || (distP2S > distP1P2)) && (distP1S < distP2S)) { // Glider not between P1 & P2 and beyond P1
                    factor = -1;
                }
                this.position = factor*distP1S/distP1P2;

                // Snap the glider point of the slider into its appropiate position
                // First, recalculate the new value of this.position
                // Second, call update(fromParent==true) to make the positioning snappier.
                if (this.snapWidth!=null && Math.abs(this._smax-this._smin)>=JXG.Math.eps) {
                    if (this.position<0.0) this.position = 0.0;
                    if (this.position>1.0) this.position = 1.0;
                    
                    var v = this.position*(this._smax-this._smin)+this._smin;
                        v = Math.round(v/this.snapWidth)*this.snapWidth;
                    this.position = (v-this._smin)/(this._smax-this._smin);
                    this.update(true);
                }
            }
            var p1Scr = this.slideObject.point1.coords.scrCoords;
            var p2Scr = this.slideObject.point2.coords.scrCoords;

            var i;
            if(this.slideObject.getSlope() == 0) {
                i = 1;
            } else {
                i = 2;
            }

            var y = this.coords.scrCoords[i];
            if(!this.slideObject.visProp['straightFirst']) {
                if(p1Scr[i] < p2Scr[i]) {
                    if(y < p1Scr[i]) {
                       this.coords = this.slideObject.point1.coords;
                       this.position = 0;
                    }
                }
                else if(p1Scr[i] > p2Scr[i]) {
                    if(y > p1Scr[i]) {
                       this.coords = this.slideObject.point1.coords;
                       this.position = 0;
                    }
                }
            }
            if(!this.slideObject.visProp['straightLast']) {
                if(p1Scr[i] < p2Scr[i]) {
                    if(y > p2Scr[i]) {
                       this.coords = this.slideObject.point2.coords;
                       this.position = 1;
                    }
                }
                else if(p1Scr[i] > p2Scr[i]) {
                    if(y < p2Scr[i]) {
                       this.coords = this.slideObject.point2.coords;
                       this.position = 1;
                    }
                }
            }  

            if(this.onPolygon) {
                var p1 = this.slideObject.point1.coords;
                var p2 = this.slideObject.point2.coords;
                if(Math.abs(this.coords.scrCoords[1]-p1.scrCoords[1])<this.board.options.precision.hasPoint && Math.abs(this.coords.scrCoords[2]-p1.scrCoords[2])<this.board.options.precision.hasPoint) {
                    var poly = this.slideObject.parentPolygon;
                    for(var i=0; i<poly.borders.length; i++) {
                        if(this.slideObject == poly.borders[i]) {
                            this.slideObject = poly.borders[(i - 1 + poly.borders.length) % poly.borders.length];
                            break;
                        }
                    }
                }
                else if(Math.abs(this.coords.scrCoords[1]-p2.scrCoords[1])<this.board.options.precision.hasPoint && Math.abs(this.coords.scrCoords[2]-p2.scrCoords[2])<this.board.options.precision.hasPoint) {
                    var poly = this.slideObject.parentPolygon;
                    for(var i=0; i<poly.borders.length; i++) {
                        if(this.slideObject == poly.borders[i]) {
                            this.slideObject = poly.borders[(i + 1 + poly.borders.length) % poly.borders.length];
                            break;                        
                        }
                    }
                }
            }
        } else if(this.slideObject.type == JXG.OBJECT_TYPE_CURVE) {
            this.updateConstraint(); // In case, the point is a constrained glider.
            this.coords  = this.board.algebra.projectPointToCurve(this, this.slideObject);
        } else if(this.slideObject.type == JXG.OBJECT_TYPE_TURTLE) {
            this.updateConstraint(); // In case, the point is a constrained glider.
            this.coords  = this.board.algebra.projectPointToTurtle(this, this.slideObject);
        }
    }
    
    /* If point is a calculated point, call updateConstraint() to calculate new coords. */
    if (this.type == JXG.OBJECT_TYPE_CAS) {
        this.updateConstraint();
    }

    this.updateTransform();
    
    //this.updateRenderer();
    this.needsUpdate = false;
    return this;
};

/**
 * Calls the renderer to update the drawing.
 * @private
 */
JXG.Point.prototype.updateRenderer = function () {
    /* Call the renderer only if point is visible. */
    if(this.visProp['visible']) {
        var wasReal = this.isReal;
        this.isReal = (isNaN(this.coords.usrCoords[1]+this.coords.usrCoords[2]))?false:true;
        this.isReal = (Math.abs(this.coords.usrCoords[0])>this.board.algebra.eps)?this.isReal:false;  //Homogeneous coords: ideal point
        if (this.isReal) {
            if (wasReal!=this.isReal) { 
                this.board.renderer.show(this); 
                if(this.hasLabel && this.label.content.visProp['visible']) this.board.renderer.show(this.label.content); 
            }
            this.board.renderer.updatePoint(this);
        } else {
            if (wasReal!=this.isReal) { 
                this.board.renderer.hide(this); 
                if(this.hasLabel && this.label.content.visProp['visible']) this.board.renderer.hide(this.label.content); 
            }
        }
    } 

    /* Update the label if visible. */
    if(this.hasLabel && this.label.content.visProp['visible'] && this.isReal) {
        //this.label.setCoordinates(this.coords);
        this.label.content.update();
        //this.board.renderer.updateLabel(this.label);
        this.board.renderer.updateText(this.label.content);
    }
    return this;
};

/**
 * Getter method for x, this is used by for CAS-points to access point coordinates.
 * @return User coordinate of point in x direction.
 * @type number
 */
JXG.Point.prototype.X = function () {
    return this.coords.usrCoords[1];
};

/**
 * Getter method for y, this is used by CAS-points to access point coordinates.
 * @return User coordinate of point in y direction.
 * @type number
 */
JXG.Point.prototype.Y = function () {
    return this.coords.usrCoords[2];
};

/**
 * Getter method for z, this is used by CAS-points to access point coordinates.
 * @return User coordinate of point in z direction.
 * @type number
 */
JXG.Point.prototype.Z = function () {
    return this.coords.usrCoords[0];
};

/**
 * New evaluation of the function term. 
 * This is required for CAS-points: Their XTerm() method is overwritten in {@link #addConstraint}
 * @return User coordinate of point in x direction.
 * @type number
 * @private
 */
JXG.Point.prototype.XEval = function () {
    return this.coords.usrCoords[1];
};

/**
 * New evaluation of the function term. 
 * This is required for CAS-points: Their YTerm() method is overwritten in {@link #addConstraint}
 * @return User coordinate of point in y direction.
 * @type number
 * @private
 */
JXG.Point.prototype.YEval = function () {
    return this.coords.usrCoords[2];
};

/**
 * New evaluation of the function term. 
 * This is required for CAS-points: Their ZTerm() method is overwritten in {@link #addConstraint}
 * @return User coordinate of point in z direction.
 * @type number
 * @private
 */
JXG.Point.prototype.ZEval = function () {
    return this.coords.usrCoords[0];
};

/**
 * Getter method for the distance to a second point, this is required for CAS-elements.
 * Here, function inlining seems to be worthwile  (for plotting).
 * @param {JXG.Point} point2 The point to which the distance shall be calculated.
 * @return Distance in user coordinate to the given point
 * @type number
 */
JXG.Point.prototype.Dist = function(point2) {
    var sum,
        c = point2.coords.usrCoords,
        ucr = this.coords.usrCoords,
        f;
        
    f = ucr[0]-c[0];
    sum = f*f;
    f = ucr[1]-c[1];
    sum += f*f;
    f = ucr[2]-c[2];
    sum += f*f;
    return Math.sqrt(sum);
    //return this.coords.distance(JXG.COORDS_BY_USER, point2.coords);
};

/**
 * Sets x and y coordinate and calls the point's update() method.
 * @param {number} method The type of coordinates used here. Possible values are {@link JXG.COORDS_BY_USER} and {@link JXG.COORDS_BY_SCREEN}.
 * @param {number} x x coordinate in screen/user units
 * @param {number} y y coordinate in screen/user units
 */
JXG.Point.prototype.setPositionDirectly = function (method, x, y) {
    var i, dx, dy, el, p,
        oldCoords = this.coords;
        
    this.coords = new JXG.Coords(method, [x,y], this.board);

    if(this.group.length != 0) {
        // Update the initial coordinates. This is needed for free points
        // that have a transformation bound to it.
        dx = this.coords.usrCoords[1]-oldCoords.usrCoords[1];
        dy = this.coords.usrCoords[2]-oldCoords.usrCoords[2];
        for (i=0;i<this.group.length;i++) {
            for (el in this.group[i].objects) {
                p = this.group[i].objects[el];
                p.initialCoords = new JXG.Coords(JXG.COORDS_BY_USER, 
                    [p.initialCoords.usrCoords[1]+dx,p.initialCoords.usrCoords[2]+dy], 
                    this.board);
            }
        }

        this.group[this.group.length-1].dX = this.coords.scrCoords[1] - oldCoords.scrCoords[1];
        this.group[this.group.length-1].dY = this.coords.scrCoords[2] - oldCoords.scrCoords[2];
        this.group[this.group.length-1].update(this);
    } else {
        // Update the initial coordinates. This is needed for free points
        // that have a transformation bound to it.
        for (i=this.transformations.length-1;i>=0;i--) {
            this.initialCoords = new JXG.Coords(method, 
                JXG.Math.matVecMult(JXG.Math.Numerics.Inverse(this.transformations[i].matrix),[1,x,y]), 
                this.board);      
        }
        this.update();
    }
    return this;
};

/**
 * TODO
 * @param {number} method The type of coordinates used here. Possible values are {@link JXG.COORDS_BY_USER} and {@link JXG.COORDS_BY_SCREEN}.
 * @param {number} x x coordinate in screen/user units
 * @param {number} y y coordinate in screen/user units
 */
JXG.Point.prototype.setPositionByTransform = function (method, x, y) {
    var oldCoords = this.coords;
    var t = this.board.create('transform',[x,y],{type:'translate'});
    if (this.transformations.length>0 && this.transformations[this.transformations.length-1].isNumericMatrix) {
        this.transformations[this.transformations.length-1].melt(t);
    } else {
        this.addTransform(this,t);
    }

    if (this.group.length != 0) {
/*
        var dCoords = new JXG.Coords(method, [x,y], this.board);
        this.group[this.group.length-1].dX = dCoords.scrCoords[1]-this.board.origin.scrCoords[1]; 
        this.group[this.group.length-1].dY = dCoords.scrCoords[2]-this.board.origin.scrCoords[2]; 
        this.group[this.group.length-1].update(this);
*/
    } else {
        this.update();
    }
    return this;
};

/**
 * Sets x and y coordinate and calls the point's update() method.
 * @param {number} method The type of coordinates used here. Possible values are {@link JXG.COORDS_BY_USER} and {@link JXG.COORDS_BY_SCREEN}.
 * @param {number} x x coordinate in screen/user units
 * @param {number} y y coordinate in screen/user units
 */
JXG.Point.prototype.setPosition = function (method, x, y) { 
    //this.setPositionByTransform(method, x, y);
    this.setPositionDirectly(method, x, y);
    return this;
};

/**
 * Convert the point to glider and update the construction.
 * @param {String,Object} glideObject The Object the point will be bound to.
 */
JXG.Point.prototype.makeGlider = function (glideObject) {
    this.slideObject = JXG.getReference(this.board, glideObject);
    this.type = JXG.OBJECT_TYPE_GLIDER;
    this.snapWidth = null;
    
    this.slideObject.addChild(this);

    if(this.slideObject.elementClass == JXG.OBJECT_CLASS_LINE) {
        this.generatePolynomial = function() {
            return this.slideObject.generatePolynomial(this);
        };
    } else if (this.slideObject.elementClass == JXG.OBJECT_CLASS_CIRCLE) {
        this.generatePolynomial = function() {
            return this.slideObject.generatePolynomial(this);
        };
    }

    //this.position = 0;
    this.needsUpdate = true;
    this.update();
    return this;
};

/**
 * Convert the point to CAS point and call update().
 * @param {array} terms [[zterm], xterm, yterm] defining terms for the z, x and y coordinate.
 * The z-coordinate is optional and it is used for homogeneaous coordinates.
 * The coordinates may be either <ul>
 *   <li>a JavaScript function,</li>
 *   <li>a string containing GEONExT syntax. This string will be converted into a JavaScript 
 *     function here,</li>
 *   <li>a number</li>
 *   <li>a pointer to a slider object. This will be converted into a call of the Value()-method 
 *     of this slider.</li>
 *   </ul>
 * @see JXG.Algebra#geonext2JS
 */
JXG.Point.prototype.addConstraint = function (terms) {
    this.type = JXG.OBJECT_TYPE_CAS;
    var elements = this.board.elementsByName;
    var newfuncs = [];
    var fs;
    
    for (var i=0;i<terms.length;i++) {
        var v = terms[i];
        if (typeof v=='string') {
            // Convert GEONExT syntax into  JavaScript syntax
            var t  = this.board.algebra.geonext2JS(v);
            newfuncs[i] = new Function('','return ' + t + ';');
        } else if (typeof v=='function') {
            newfuncs[i] = v;
        } else if (typeof v=='number') {
            newfuncs[i] = function(z){ return function() { return z; }; }(v);
        } else if (typeof v == 'object' && typeof v.Value == 'function') {    // Slider
            newfuncs[i] = (function(a) { return function() { return a.Value(); };})(v);
        }
    }
    if (terms.length==1) { // Intersection function
        this.updateConstraint = function() { 
                var c = newfuncs[0](); 
                if (JXG.isArray(c)) {      // Array
                    this.coords.setCoordinates(JXG.COORDS_BY_USER,c);
                } else {                   // Coords object
                    this.coords = c;
                }
            };
        // if (!this.board.isSuspendedUpdate) { this.update(); }
        // return this;
    } else if (terms.length==2) { // Euclidean coordinates
        this.XEval = newfuncs[0];
        this.YEval = newfuncs[1];
        fs = 'this.coords.setCoordinates(JXG.COORDS_BY_USER,[this.XEval(),this.YEval()]);';
        this.updateConstraint = new Function('',fs);
    } else { // Homogeneous coordinates
        this.ZEval = newfuncs[0];
        this.XEval = newfuncs[1];
        this.YEval = newfuncs[2];
        fs = 'this.coords.setCoordinates(JXG.COORDS_BY_USER,[this.ZEval(),this.XEval(),this.YEval()]);';
        this.updateConstraint = new Function('',fs);
    }

    if (!this.board.isSuspendedUpdate) { this.update(); }
    return this;
};

/**
 * TODO
 */
JXG.Point.prototype.updateTransform = function () {
    if (this.transformations.length==0 || this.baseElement==null) {
        return;
    }
    var c, i;

    if (this===this.baseElement) {      // case of bindTo
        c = this.transformations[0].apply(this.baseElement,'self');
    } else {                           // case of board.create('point',[baseElement,transform]);
        c = this.transformations[0].apply(this.baseElement);
    }
    this.coords.setCoordinates(JXG.COORDS_BY_USER,c);
    for (i=1;i<this.transformations.length;i++) {
        this.coords.setCoordinates(JXG.COORDS_BY_USER,this.transformations[i].apply(this));
    }
    return this;
};

/**
 * TODO
 * @param el TODO
 * @param transform TODO
 */
JXG.Point.prototype.addTransform = function (el, transform) {
    var list, i, len;
    if (this.transformations.length==0) { // There is only one baseElement possible
        this.baseElement = el;
    }
    if (JXG.isArray(transform)) {
        list = transform;
    } else {
        list = [transform];
    }
    len = list.length;
    for (i=0;i<len;i++) {
        this.transformations.push(list[i]);
    }
    return this;
};

/**
 * Animate the point. 
 * @param {number} direction The direction the glider is animated. Can be +1 or -1.
 * @param {number} stepCount The number of steps.
 * @name Glider#startAnimation
 * @see Glider#stopAnimation
 * @function
 */
JXG.Point.prototype.startAnimation = function(direction, stepCount) {
    if((this.type == JXG.OBJECT_TYPE_GLIDER) && (typeof this.intervalCode == 'undefined')) {
        this.intervalCode = window.setInterval('JXG.JSXGraph.boards[\'' + this.board.id + '\'].objects[\'' + this.id + '\']._anim(' + direction + ', ' + stepCount + ')', 250);
        if(typeof this.intervalCount == 'undefined')
            this.intervalCount = 0;
    }
    return this;
};

/**
 * Stop animation.
 * @name Glider#stopAnimation
 * @see Glider#startAnimation
 * @function
 */
JXG.Point.prototype.stopAnimation = function() {
    if(typeof this.intervalCode != 'undefined') {
        window.clearInterval(this.intervalCode);
        delete(this.intervalCode);
    }
    return this;
};

/**
 * Starts an animated point movement towards the given coordinates <tt>where</tt>. The animation is done after <tt>time</tt> milliseconds.
 * @param {Array} where Array containing the x and y coordinate of the target location.
 * @param {int} time Number of milliseconds the animation should last.
 * If the second parameter is not given or is equal to 0, setPosition() is called, see #setPosition.
 * @see #animate
 */
JXG.Point.prototype.moveTo = function(where, time) {
    if (typeof time == 'undefined' || time == 0) {
        this.setPosition(JXG.COORDS_BY_USER, where[0], where[1]);
        //this.prepareUpdate().update().updateRenderer();
        this.board.update(this);
        return this;
    }
	var delay = 35,
	    steps = Math.ceil(time/(delay * 1.0)),
		coords = new Array(steps+1),
		X = this.coords.usrCoords[1],
		Y = this.coords.usrCoords[2],
		dX = (where[0] - X),
		dY = (where[1] - Y),
	    i;
    
    if(Math.abs(dX) < JXG.Math.eps && Math.abs(dY) < JXG.Math.eps)
        return this;
	
	for(i=steps; i>=0; i--) {
		coords[steps-i] = [X + dX * Math.sin((i/(steps*1.0))*Math.PI/2.), Y+ dY * Math.sin((i/(steps*1.0))*Math.PI/2.)];
	}
	this.animationPath = coords;
	this.board.animationObjects[this.id] = this;
	if(typeof this.board.animationIntervalCode == 'undefined') {
		this.board.animationIntervalCode = window.setInterval('JXG.JSXGraph.boards[\'' + this.board.id + '\'].animate();', delay);
	}
    return this;
};

/**
 * Starts an animated point movement towards the given coordinates <tt>where</tt>. After arriving at <tt>where</tt> the point moves back to where it started.
 * The animation is done after <tt>time</tt> milliseconds.
 * @param {Array} where Array containing the x and y coordinate of the target location.
 * @param {int} time Number of milliseconds the animation should last.
 * @param {int} repeat Optional: How often the animation should be repeated. The time value is then taken for one repeat.
 * @see #animate
 */
JXG.Point.prototype.visit = function(where, time, repeat) {
    if(arguments.length == 2)
        repeat = 1;

    var delay = 35,
        steps = Math.ceil(time/(delay * 1.0)),
        coords = new Array(repeat*(steps+1)),
        X = this.coords.usrCoords[1],
        Y = this.coords.usrCoords[2],
        dX = (where[0] - X),
        dY = (where[1] - Y),
        i, j;
    
    for(j=0; j<repeat; j++) {
        for(i=steps; i>=0; i--) {
            coords[j*(steps+1) + steps-i] = [X + dX * Math.pow(Math.sin((i/(steps*1.0))*Math.PI), 2.), Y+ dY * Math.pow(Math.sin((i/(steps*1.0))*Math.PI), 2.)];
        }
    }
    this.animationPath = coords;
    this.board.animationObjects[this.id] = this;
    if(typeof this.board.animationIntervalCode == 'undefined') {
        this.board.animationIntervalCode = window.setInterval('JXG.JSXGraph.boards[\'' + this.board.id + '\'].animate();', delay);
    }
    return this;
};

/**
 * Animates a glider. Is called by the browser after startAnimation is called.
 * @param {number} direction The direction the glider is animated.
 * @param {number} stepCount The number of steps.
 * @see #startAnimation
 * @see #stopAnimation
 * @private
 */
JXG.Point.prototype._anim = function(direction, stepCount) {
    var distance, slope, dX, dY, alpha, startPoint,
        factor = 1, newX, radius;
    
    this.intervalCount++;
    if(this.intervalCount > stepCount)
        this.intervalCount = 0;
    
    if(this.slideObject.type == JXG.OBJECT_TYPE_LINE) {
        distance = this.slideObject.point1.coords.distance(JXG.COORDS_BY_SCREEN, this.slideObject.point2.coords);
        slope = this.slideObject.getSlope();
        if(slope != 'INF') {
            alpha = Math.atan(slope);
            dX = Math.round((this.intervalCount/stepCount) * distance*Math.cos(alpha));
            dY = Math.round((this.intervalCount/stepCount) * distance*Math.sin(alpha));
        } else {
            dX = 0;
            dY = Math.round((this.intervalCount/stepCount) * distance);
        }
        
        if(direction < 0) {
            startPoint = this.slideObject.point2;
            if(this.slideObject.point2.coords.scrCoords[1] - this.slideObject.point1.coords.scrCoords[1] > 0)
                factor = -1;
            else if(this.slideObject.point2.coords.scrCoords[1] - this.slideObject.point1.coords.scrCoords[1] == 0) {
                if(this.slideObject.point2.coords.scrCoords[2] - this.slideObject.point1.coords.scrCoords[2] > 0)
                    factor = -1;
            }
        } else {
            startPoint = this.slideObject.point1;
            if(this.slideObject.point1.coords.scrCoords[1] - this.slideObject.point2.coords.scrCoords[1] > 0)
                factor = -1;
            else if(this.slideObject.point1.coords.scrCoords[1] - this.slideObject.point2.coords.scrCoords[1] == 0) {
                if(this.slideObject.point1.coords.scrCoords[2] - this.slideObject.point2.coords.scrCoords[2] > 0)
                    factor = -1;
            }
        }
        
        this.coords.setCoordinates(JXG.COORDS_BY_SCREEN, [startPoint.coords.scrCoords[1] + factor*dX, startPoint.coords.scrCoords[2] + factor*dY]);
    } else if(this.slideObject.type == JXG.OBJECT_TYPE_CURVE) {
        if(direction > 0) {
            newX = Math.round(this.intervalCount/stepCount * this.board.canvasWidth);
        } else {
            newX = Math.round((stepCount - this.intervalCount)/stepCount * this.board.canvasWidth);
        }
  
        this.coords.setCoordinates(JXG.COORDS_BY_SCREEN, [newX, 0]);
        this.coords = this.board.algebra.projectPointToCurve(this, this.slideObject);
    } else if(this.slideObject.type == JXG.OBJECT_TYPE_CIRCLE) {
        if(direction < 0) {
            alpha = this.intervalCount/stepCount * 2*Math.PI;
        } else {
            alpha = (stepCount - this.intervalCount)/stepCount * 2*Math.PI;
        }

        radius = this.slideObject.Radius();

        this.coords.setCoordinates(JXG.COORDS_BY_USER, [this.slideObject.midpoint.coords.usrCoords[1] + radius*Math.cos(alpha), this.slideObject.midpoint.coords.usrCoords[2] + radius*Math.sin(alpha)]);
    }
    
    this.board.update(this);
    return this;
};

/**
 * Set the style of a point.
 * @param {int} i Integer to determine the style. See {@link JXG.GeometryElement#style} for a list of available styles.
 * @see JXG.GeometryElement#style
 * @private
 * @deprecated
 */
JXG.Point.prototype.setStyle = function(i) {
    if(i == 0 || i == 1 || i == 2) { // x
        this.visProp['face'] = 'cross';
        if(i == 0) {
            this.visProp['size'] = 2;
        }
        else if(i == 1) {
            this.visProp['size'] = 3;
        }
        else {
            this.visProp['size'] = 4;
        }        
    }
    else if(i == 3 || i == 4 || i == 5 || i == 6) { // circle
        this.visProp['face'] = 'circle';
        if(i == 3) {
            this.visProp['size'] = 1;
        }
        else if(i == 4) {
            this.visProp['size'] = 2;
        }
        else if(i == 5) {
            this.visProp['size'] = 3;
        }        
        else {
            this.visProp['size'] = 4;
        }            
    }
    else if(i == 7 || i == 8 || i == 9) { // rectangle
        this.visProp['face'] = 'square';
        if(i == 7) {
            this.visProp['size'] = 2;
        }
        else if(i == 8) {
            this.visProp['size'] = 3;
        }
        else {
            this.visProp['size'] = 4;
        }  
    }
    else if(i == 10 || i == 11 || i == 12) { // +
        this.visProp['face'] = 'plus';
        if(i == 10) {
            this.visProp['size'] = 2;
        }
        else if(i == 11) {
            this.visProp['size'] = 3;
        }
        else {
            this.visProp['size'] = 4;
        }  
    }    
    
    this.board.renderer.changePointStyle(this);
    return this;
};

/**
 * Set the face of a point.
 * @param {string} s String which determines the face of the point. See {@link JXG.GeometryElement#face} for a list of available faces.
 * @see JXG.GeometryElement#face
 * @private
 */
JXG.Point.prototype.setFace = function(s) {
    s = s.toLowerCase();
    if(s == 'cross' || s == 'x' || s == 'plus' || s == '+' || s == 'circle' || s == 'o' || s == 'square' || s == '[]' 
       || s == 'diamond' || s == '<>' || s == 'triangleup' || s == 'a' || s == 'triangledown' || s == 'v' || 
       s == 'triangleleft' || s == '<' || s == 'triangleright' || s == '>') {
        this.visProp['face'] = s;
    }
    else {
        this.visProp['face'] = 'circle';
    }
    this.board.renderer.changePointStyle(this);
    return this;
};

/**
 * Remove the point from the drawing.
 */
JXG.Point.prototype.remove = function() {    
    if (this.hasLabel) {
        this.board.renderer.remove(document.getElementById(this.label.content.id));
    }
    this.board.renderer.remove(document.getElementById(this.id));
};

/**
 * TODO
 * @return TODO
 * @type JXG.Coords
 * @private
 */
JXG.Point.prototype.getTextAnchor = function() {
    return this.coords;
};

/**
 * TODO
 * @return TODO
 * @type JXG.Coords
 * @private
 */
JXG.Point.prototype.getLabelAnchor = function() {
    return this.coords;
};

/**
 * Copy the element to the background.
 * @param addToTrace If true the clone will be added to trace control and can be removed using {@link JXG.GeometryElement#clearTrace}.
 * Currently not used, and always true.
 */
JXG.Point.prototype.cloneToBackground = function(/** boolean */ addToTrace) {
    var copy = {};
    copy.id = this.id + 'T' + this.numTraces;
    this.numTraces++;
    copy.coords = this.coords;
    copy.visProp = this.visProp;
    copy.elementClass = JXG.OBJECT_CLASS_POINT;
    JXG.clearVisPropOld(copy);
    
    this.board.renderer.drawPoint(copy);

    this.traces[copy.id] = document.getElementById(copy.id);

    delete copy;
/*   
    this.board.renderer.cloneSubTree(this);
*/    
    return this;
};

/* old description of the following createPoint method.
 * There are several methods to construct a point.
 * The input parameter "parentArr" determines the point:
 * - 2 numbers: affine (Euclidean) coordinates of a free point
 * - 2 numbers and atts['slideObject'] : Glider with initial Euclidean coordinates
 * - 2 Strings or (1 String and 1 Number): constrained point
 * - 1 function: intersection of objects, this is just a constrained point too
 * - 1 transformation object: clone of a base point transformed by the given Transformation
 * - 3 numbers: homogeneous coordinates of a free point
 */

/**
 * @class This element is used to provide a constructor for a general point. A free point is created if the given parent elements are all numbers
 * and the property fixed is not set or set to false. If one or more parent elements is not a number but a string containing a GEONE<sub>x</sub>T
 * constraint or a function the point will be considered as constrained). That means that the user won't be able to change the point's
 * position directly.
 * @pseudo
 * @description
 * @name Point
 * @augments JXG.Point
 * @constructor
 * @type JXG.Point
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {number,string,function_number,string,function_number,string,function} z_,x,y Parent elements can be two or three elements of type number, a string containing a GEONE<sub>x</sub>T
 * constraint, or a function which takes no parameter and returns a number. Every parent element determines one coordinate. If a coordinate is
 * given by a number, the number determines the initial position of a free point. If given by a string or a function that coordinate will be constrained
 * that means the user won't be able to change the point's position directly by mouse because it will be calculated automatically depending on the string
 * or the function's return value. If two parent elements are given the coordinates will be interpreted as 2D affine euclidean coordinates, if three such
 * parent elements are given they will be interpreted as homogeneous coordinates.
 * @param {JXG.Point_JXG.Transformation} Point,Transformation A point can also be created providing a transformation. The resulting point is a clone of the base
 * point transformed by the given Transformation. {@see JXG.Transformation}.
 * @example
 * // Create a free point using affine euclidean coordinates 
 * var p1 = board.create('point', [3.5, 2.0]);
 * </pre><div id="672f1764-7dfa-4abc-a2c6-81fbbf83e44b" style="width: 200px; height: 200px;"></div>
 * <script type="text/javascript">
 *   var board = JXG.JSXGraph.initBoard('672f1764-7dfa-4abc-a2c6-81fbbf83e44b', {boundingbox: [-1, 5, 5, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var p1 = board.create('point', [3.5, 2.0]);
 * </script><pre>
 * @example
 * // Create a constrained point using anonymous function 
 * var p2 = board.create('point', [3.5, function () { return p1.X(); }]);
 * </pre><div id="4fd4410c-3383-4e80-b1bb-961f5eeef224" style="width: 200px; height: 200px;"></div>
 * <script type="text/javascript">
 *   var fpex1_board = JXG.JSXGraph.initBoard('4fd4410c-3383-4e80-b1bb-961f5eeef224', {boundingbox: [-1, 5, 5, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var fpex1_p1 = fpex1_board.create('point', [3.5, 2.0]);
 *   var fpex1_p2 = fpex1_board.create('point', [3.5, function () { return fpex1_p1.X(); }]);
 * </script><pre>
 * @example
 * // Create a point using transformations 
 * var trans = board.create('transform', [2, 0.5], {type:'scale'});
 * var p3 = board.create('point', [p2, trans]);
 * </pre><div id="630afdf3-0a64-46e0-8a44-f51bd197bb8d" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var fpex2_board = JXG.JSXGraph.initBoard('630afdf3-0a64-46e0-8a44-f51bd197bb8d', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var fpex2_trans = fpex2_board.create('transform', [2, 0.5], {type:'scale'});
 *   var fpex2_p2 = fpex2_board.create('point', [3.5, 2.0]);
 *   var fpex2_p3 = fpex2_board.create('point', [fpex2_p2, fpex2_trans]);
 * </script><pre>
 */
JXG.createPoint = function(/** JXG.Board */ board, /** array */ parents, /** object */ atts) {
    var el, isConstrained = false, i, show;
    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'point','withLabel'), layer:null});
    show = (typeof atts['visible']=='undefined') || JXG.str2Bool(atts['visible']);
    
    for (i=0;i<parents.length;i++) {
        if (typeof parents[i]=='function' || typeof parents[i]=='string') {
            isConstrained = true;
        }
    }
    if (!isConstrained) {
        if ( (JXG.isNumber(parents[0])) && (JXG.isNumber(parents[1])) ) {
            el = new JXG.Point(board, parents, atts['id'], atts['name'], show, atts['withLabel'], atts['layer']);
            if ( atts["slideObject"] != null ) {
                el.makeGlider(atts["slideObject"]);
            } else {
                el.baseElement = el; // Free point
            }
        } else if ( (typeof parents[0]=='object') && (typeof parents[1]=='object') ) { // Transformation
            el = new JXG.Point(board, [0,0], atts['id'], atts['name'], show, atts['withLabel'], atts['layer']);   
            el.addTransform(parents[0],parents[1]);
        }
        else {// Failure
            throw new Error("JSXGraph: Can't create point with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
        }
    } else {
        el = new JXG.Point(board, [0,0], atts['id'], atts['name'], show, atts['withLabel'], atts['layer']);
        el.addConstraint(parents);
    }
    return el;
};

/**
 * @class This element is used to provide a constructor for a glider point. 
 * @pseudo
 * @description A glider is a point which lives on another geometric element like a line, circle, curve, turtle.
 * @name Glider
 * @augments JXG.Point
 * @constructor
 * @type JXG.Point
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {number_number_number_JXG.GeometryElement} z_,x_,y_,GlideObject Parent elements can be two or three elements of type number and the object the glider lives on.
 * The coordinates are completely optional. If not given the origin is used. If you provide two numbers for coordinates they will be interpreted as affine euclidean
 * coordinates, otherwise they will be interpreted as homogeneous coordinates. In any case the point will be projected on the glide object.
 * @example
 * // Create a glider with user defined coordinates. If the coordinates are not on
 * // the circle (like in this case) the point will be projected onto the circle.
 * var p1 = board.create('point', [2.0, 2.0]);
 * var c1 = board.create('circle', [p1, 2.0]);
 * var p2 = board.create('glider', [2.0, 1.5, c1]);
 * </pre><div id="4f65f32f-e50a-4b50-9b7c-f6ec41652930" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var gpex1_board = JXG.JSXGraph.initBoard('4f65f32f-e50a-4b50-9b7c-f6ec41652930', {boundingbox: [-1, 5, 5, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var gpex1_p1 = gpex1_board.create('point', [2.0, 2.0]);
 *   var gpex1_c1 = gpex1_board.create('circle', [gpex1_p1, 2.0]);
 *   var gpex1_p2 = gpex1_board.create('glider', [2.0, 1.5, gpex1_c1]);
 * </script><pre>
 * @example
 * // Create a glider with default coordinates (1,0,0). Same premises as above.
 * var p1 = board.create('point', [2.0, 2.0]);
 * var c1 = board.create('circle', [p1, 2.0]);
 * var p2 = board.create('glider', [c1]);
 * </pre><div id="4de7f181-631a-44b1-a12f-bc4d995609e8" style="width: 200px; height: 200px;"></div>
 * <script type="text/javascript">
 *   var gpex2_board = JXG.JSXGraph.initBoard('4de7f181-631a-44b1-a12f-bc4d995609e8', {boundingbox: [-1, 5, 5, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var gpex2_p1 = gpex2_board.create('point', [2.0, 2.0]);
 *   var gpex2_c1 = gpex2_board.create('circle', [gpex2_p1, 2.0]);
 *   var gpex2_p2 = gpex2_board.create('glider', [gpex2_c1]);
 * </script><pre>
 */
JXG.createGlider = function(board, parents, atts) {
    var el, show;
    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'point','withLabel'), layer:null});
    show = (typeof atts['visible']=='undefined') || JXG.str2Bool(atts['visible']);
    
    if (parents.length==1) {
      el = new JXG.Point(board, [0,0], atts['id'], atts['name'], show, atts['withLabel']);
    } else {
      el = board.create('point',parents.slice(0,-1), atts);
    }
    el.makeGlider(parents[parents.length-1]);
    return el;
};

/**
 * @class This element is used to provide a constructor for an intersection point. 
 * @pseudo
 * @description An intersection point is a point which lives on two Lines or Circles or one Line and one Circle at the same time, i.e.
 * an intersection point of the two elements.
 * @name Intersection
 * @augments JXG.Point
 * @constructor
 * @type JXG.Point
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Line,JXG.Circle_JXG.Line,JXG.Circle_number} el1,el2,i The result will be a intersection point on el1 and el2. i determines the
 * intersection point if two points are available: <ul>
 *   <li>i==0: use the positive square root,</li> 
 *   <li>i==1: use the negative square root.</li></ul>
 * @example
 * // Create an intersection point of circle and line
 * var p1 = board.create('point', [2.0, 2.0]);
 * var c1 = board.create('circle', [p1, 2.0]);
 * 
 * var p2 = board.create('point', [2.0, 2.0]);
 * var p3 = board.create('point', [2.0, 2.0]);
 * var l1 = board.create('line', [p2, p3]);
 * 
 * var i = board.create('intersection', [c1, l1, 0]);
 * </pre><div id="e5b0e190-5200-4bc3-b995-b6cc53dc5dc0" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var ipex1_board = JXG.JSXGraph.initBoard('e5b0e190-5200-4bc3-b995-b6cc53dc5dc0', {boundingbox: [-1, 7, 7, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var ipex1_p1 = ipex1_board.create('point', [4.0, 4.0]);
 *   var ipex1_c1 = ipex1_board.create('circle', [ipex1_p1, 2.0]);
 *   var ipex1_p2 = ipex1_board.create('point', [1.0, 1.0]);
 *   var ipex1_p3 = ipex1_board.create('point', [5.0, 3.0]);
 *   var ipex1_l1 = ipex1_board.create('line', [ipex1_p2, ipex1_p3]);
 *   var ipex1_i = ipex1_board.create('intersection', [ipex1_c1, ipex1_l1, 0]);
 * </script><pre>
 */
JXG.createIntersectionPoint = function(board, parents, attributes) {
    var el;
    if (parents.length>=3) {
        if(parents.length == 3)
            parents.push(null);
        el = board.create('point', [board.intersection(parents[0], parents[1], parents[2], parents[3])], attributes);
    }

    parents[0].addChild(el);
    parents[1].addChild(el);

    el.generatePolynomial = function () {
        var poly1 = parents[0].generatePolynomial(el);
        var poly2 = parents[1].generatePolynomial(el);

        if((poly1.length == 0) || (poly2.length == 0))
            return [];
        else
            return [poly1[0], poly2[0]];
    };
    
    return el;
};

/**
 * @class This element is used to provide a constructor for the "other" intersection point.
 * @pseudo
 * @description An intersection point is a point which lives on two Lines or Circles or one Line and one Circle at the same time, i.e.
 * an intersection point of the two elements. Additionally, one intersection point is provided. The function returns the other intersection point.
 * @name OtherIntersection
 * @augments JXG.Point
 * @constructor
 * @type JXG.Point
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Line,JXG.Circle_JXG.Line,JXG.Circle_JXG.Point} el1,el2,p The result will be a intersection point on el1 and el2. i determines the
 * intersection point different from p: 
 * @example
 * // Create an intersection point of circle and line
 * var p1 = board.create('point', [2.0, 2.0]);
 * var c1 = board.create('circle', [p1, 2.0]);
 * 
 * var p2 = board.create('point', [2.0, 2.0]);
 * var p3 = board.create('point', [2.0, 2.0]);
 * var l1 = board.create('line', [p2, p3]);
 * 
 * var i = board.create('intersection', [c1, l1, 0]);
 * var j = board.create('otherintersection', [c1, l1, i]);
 * </pre><div id="45e25f12-a1de-4257-a466-27a2ae73614c" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var ipex2_board = JXG.JSXGraph.initBoard('45e25f12-a1de-4257-a466-27a2ae73614c', {boundingbox: [-1, 7, 7, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var ipex2_p1 = ipex1_board.create('point', [4.0, 4.0]);
 *   var ipex2_c1 = ipex1_board.create('circle', [ipex2_p1, 2.0]);
 *   var ipex2_p2 = ipex1_board.create('point', [1.0, 1.0]);
 *   var ipex2_p3 = ipex1_board.create('point', [5.0, 3.0]);
 *   var ipex2_l1 = ipex1_board.create('line', [ipex2_p2, ipex2_p3]);
 *   var ipex2_i = ipex1_board.create('intersection', [ipex2_c1, ipex2_l1, 0]);
 *   var ipex2_j = ipex1_board.create('intersection', [ipex2_c1, ipex2_l1, ipex2_i]);
 * </script><pre>
 */
JXG.createOtherIntersectionPoint = function(board, parents, attributes) {
    var el;
    if (parents.length!=3 || 
        !JXG.isPoint(parents[2]) ||
        (parents[0].elementClass != JXG.OBJECT_CLASS_LINE && parents[0].elementClass != JXG.OBJECT_CLASS_CIRCLE) ||
        (parents[1].elementClass != JXG.OBJECT_CLASS_LINE && parents[1].elementClass != JXG.OBJECT_CLASS_CIRCLE) ) {
        // Failure
        throw new Error("JSXGraph: Can't create 'other intersection point' with parent types '" + (typeof parents[0]) + "',  '" + (typeof parents[1])+ "'and  '" + (typeof parents[2]) + "'.");
    }
    else {
        el = board.create('point', [board.otherIntersection(parents[0], parents[1], parents[2])], attributes);
    }
    
    parents[0].addChild(el);
    parents[1].addChild(el);

    el.generatePolynomial = function () {
        var poly1 = parents[0].generatePolynomial(el);
        var poly2 = parents[1].generatePolynomial(el);

        if((poly1.length == 0) || (poly2.length == 0))
            return [];
        else
            return [poly1[0], poly2[0]];
    };
    
    return el;
};


JXG.JSXGraph.registerElement('point', {
    icon:           'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAC0AAAAqCAIAAACofUV1AAAAu0lEQVR42u2Y2w7EIAhEYeL//zL7ttntRakyhjTw2MhwYDSaqplJgoDkiAAOFV0XaSGFD19MjMjh7/u70g8E6vBV1JmIQK2VHhp7A/5KdWxKf24Dh+HRxDaIvnJiX3jD6OhjM8RdlRfdc/Ece0y5rFW+FEdxFMerOCbe+9NxqFW+9Dn2WHOuAs8iNkT6/cEbyZ0yniYwIBL50ob4IY/F4XThkVj0yJOOQK2VHtpEW0OnuP+lqEep7rn/+AD75zNf8mTQTQAAAABJRU5ErkJggg%3D%3D',
    label:          'Free point',
    alttext:        'Constructs a free point',
    category:       'basic/points',
    description:    'Click on the board to place a free point or enter a pair of coordinates in the textbox.',
    showCoordsBox:  true,
    showInputbox:   false,
    checkInput:     function (draft, input) {
                       if(draft && input[input.length-1].usrCoords)
                           return true;

                       if(!draft && input.length == 1) {
                           return board.create('point', input[0].usrCoords.slice(1));
                       }

                       return false;
                    },
    creator:        JXG.createPoint
});
JXG.JSXGraph.registerElement('glider', JXG.createGlider);
JXG.JSXGraph.registerElement('intersection', JXG.createIntersectionPoint);
JXG.JSXGraph.registerElement('otherintersection', JXG.createOtherIntersectionPoint);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview The geometry object Line is defined in this file. Line stores all
 * style and functional properties that are required to draw and move a line on
 * a board.
 */

/**
 * The Line class is a basic class for all kind of line objects, e.g. line, arrow, and axis. It is usually defined by two points and can
 * be intersected with some other geometry elements.
 * @class Creates a new basic line object. Do not use this constructor to create a line. Use {@link JXG.Board#create} with
 * type {@link Line}, {@link Arrow}, or {@link Axis} instead.
 * @constructor
 * @augments JXG.GeometryElement
 * @param {String,JXG.Board} board The board the new line is drawn on.
 * @param {Point} p1 Startpoint of the line.
 * @param {Point} p2 Endpoint of the line.
 * @param {String} id Unique identifier for this object. If null or an empty string is given,
 * an unique id will be generated by Board
 * @param {String} name Not necessarily unique name. If null or an
 * empty string is given, an unique name will be generated.
 * @param {boolean} withLabel construct label, yes/no
 * @param {integer} layer display layer [0-9]
 * @see JXG.Board#generateName
 */
JXG.Line = function (board, p1, p2, id, name, withLabel, layer) {
    /* Call the constructor of GeometryElement */
    this.constructor();

   /**
     * Sets type of GeometryElement, value is OBJECT_TYPE_LINE.
     * @constant
     * @type int
     * @default JXG#OBJECT_TYPE_LINE
     * @private
     */
    this.type = JXG.OBJECT_TYPE_LINE;

    /**
     * Class of element, value is OBJECT_CLASS_LINE;
     * @type int
     * @constant
     * @default JXG#OBJECT_CLASS_LINE
     * @private
     */
    this.elementClass = JXG.OBJECT_CLASS_LINE;

    this.init(board, id, name);

    /**
     * Set the display layer.
     */
    if (layer == null) layer = board.options.layer['line'];
    this.layer = layer;
    
    /**
     * Startpoint of the line. You really should not set this field directly as it may break JSXGraph's
     * udpate system so your construction won't be updated properly.
     * @type JXG.Point
     */
    this.point1 = JXG.getReference(this.board, p1);

    /**
     * Endpoint of the line. Just like {@link #point1} you shouldn't write this field directly.
     * @type JXG.Point
     */
    this.point2 = JXG.getReference(this.board, p2);

    /**
     * An image bound to this line.
     * @type JXG.Image
     * @default null
     * @private
     */
    this.image = null;

    /**
     * TODO description
     * @type array
     * @default [[1,0,0],[0,1,0],[0,0,1]]
     * @private
     */
    this.imageTransformMatrix = [[1,0,0],[0,1,0],[0,0,1]];

    /**
     * This is just for the hasPoint() method.
     * @type int
     * @private
     */
    //this.r = this.board.options.precision.hasPoint;

    /* Wurde in GeometryElement schon dokumentiert. */
    this.visProp['fillColor'] = this.board.options.line.fillColor;
    this.visProp['highlightFillColor'] = this.board.options.line.highlightFillColor;
    this.visProp['strokeColor'] = this.board.options.line.strokeColor;
    this.visProp['highlightStrokeColor'] = this.board.options.line.highlightStrokeColor;

    /**
     * Determines if a line is drawn beyond {@link #point1}.
     * @name JXG.Line#straightFirst
     * @type boolean
     * @default JXG.Options.line#straightFirst
     * @field
     * @see JXG.Line#straightLast
     */
    this.visProp['straightFirst'] = this.board.options.line.straightFirst;

    /**
     * Determines if a line is drawn beyond {@link #point2}.
     * @name JXG.Line#straightLast
     * @type boolean
     * @default JXG.Options.line#straightLast
     * @field
     * @see JXG.Line#straightFirst
     */
    this.visProp['straightLast'] = this.board.options.line.straightLast;

    /* Already documented in JXG.GeometryElement */
    this.visProp['visible'] = true;

    /**
     * Determines if a line has an arrow at {@link #point1}.
     * @name JXG.Line#firstArrow
     * @type boolean
     * @default JXG.Options.line#firstArrow
     * @field
     * @see JXG.Line#lastArrow
     */
    this.visProp['firstArrow'] = this.board.options.line.firstArrow;

    /**
     * Determines if a line has an arrow at {@link #point1}.
     * @name JXG.Line#lastArrow
     * @type boolean
     * @default JXG.Options.line#lastArrow
     * @field
     * @see JXG.Line#firstArrow
     */
    this.visProp['lastArrow'] = this.board.options.line.lastArrow;

    /**
     * Array of ticks storing all the ticks on this line. Do not set this field directly and use
     * {@link #addTicks} and {@link #removeTicks} to add and remove ticks to and from the line.
     * @type array
     * @see JXG.Ticks
     */
    this.ticks = [];

    /**
     * Reference of the ticks created automatically when constructing an axis.
     * @type JXG.Ticks
     * @see JXG.Ticks
     */
    this.defaultTicks = null;

    /**
    * If the line is the border of a polygon, the polygon object is stored, otherwise null.
    * @type JXG.Polygon
    * @default null
    * @private
    */
    this.parentPolygon = null;

    // create Label
    this.createLabel(withLabel);

    /* Register line at board */
    this.id = this.board.addLine(this);

    /* Add arrow as child to defining points */
    this.point1.addChild(this);
    this.point2.addChild(this);
    this.update();
};

JXG.Line.prototype = new JXG.GeometryElement;

/**
 * Checks whether (x,y) is near the line.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @return {boolean} True if (x,y) is near the line, False otherwise.
 */
 JXG.Line.prototype.hasPoint = function (x, y) {
    // Compute the stdform of the line in screen coordinates.
    var c = [], s,
        v = [1,x,y],
        vnew = [],
        mu, i, coords, p1Scr, p2Scr, distP1P, distP2P, distP1P2;

    c[0] = this.stdform[0] -
                this.stdform[1]*this.board.origin.scrCoords[1]/this.board.stretchX+
                this.stdform[2]*this.board.origin.scrCoords[2]/this.board.stretchY;
    c[1] = this.stdform[1]/this.board.stretchX;
    c[2] = this.stdform[2]/(-this.board.stretchY);

    // Project the point orthogonally onto the line 
    var vnew = [0,c[1],c[2]];
    vnew = JXG.Math.crossProduct(vnew,v); // Orthogonal line to c through v
    vnew = JXG.Math.crossProduct(vnew,c); // Intersect orthogonal line with line

    // Normalize the projected point
    vnew[1] /= vnew[0];
    vnew[2] /= vnew[0];
    vnew[0] = 1.0;
    
    // The point is too far away from the line
    // dist(v,vnew)^2 projective
    //if (this.board.algebra.distance(v,vnew)>this.board.options.precision.hasPoint) {
    s = (v[0]-vnew[0])*(v[0]-vnew[0])+(v[1]-vnew[1])*(v[1]-vnew[1])+(v[2]-vnew[2])*(v[2]-vnew[2]);
    if (isNaN(s) || s>this.board.options.precision.hasPoint*this.board.options.precision.hasPoint) {
        return false;
    }

    if(this.visProp['straightFirst'] && this.visProp['straightLast']) {
        return true;
    } else { // If the line is a ray or segment we have to check if the projected point is "inside" P1 and P2.
/*
        coords = new JXG.Coords(JXG.COORDS_BY_SCREEN, [vnew[1],vnew[2]], this.board);
        p1Scr = this.point1.coords.scrCoords;
        p2Scr = this.point2.coords.scrCoords;
        distP1P = coords.distance(JXG.COORDS_BY_SCREEN, this.point1.coords);
        distP2P = coords.distance(JXG.COORDS_BY_SCREEN, this.point2.coords);
        distP1P2 = this.point1.coords.distance(JXG.COORDS_BY_SCREEN, this.point2.coords);
*/
        p1Scr = this.point1.coords.scrCoords;
        p2Scr = this.point2.coords.scrCoords;
        distP1P2 = (p2Scr[1]-p1Scr[1])*(p2Scr[1]-p1Scr[1])+(p2Scr[2]-p1Scr[2])*(p2Scr[2]-p1Scr[2]);  // dist(p1,p2)^2 affine
        distP1P = (vnew[1]-p1Scr[1])*(vnew[1]-p1Scr[1])+(vnew[2]-p1Scr[2])*(vnew[2]-p1Scr[2]);       // dist(vnew,p1)^2 affine
        distP2P = (vnew[1]-p2Scr[1])*(vnew[1]-p2Scr[1])+(vnew[2]-p2Scr[2])*(vnew[2]-p2Scr[2]);       // dist(vnew,p2)^2 affine

        if((distP1P > distP1P2) || (distP2P > distP1P2)) { // Check if P(x|y) is not between  P1 and P2
            if(distP1P < distP2P) { // P liegt auf der Seite von P1
                if(!this.visProp['straightFirst']) {
                    return false;
                }
            } else { // P liegt auf der Seite von P2
                if(!this.visProp['straightLast']) {
                    return false;
                }
            }
        }
        return true;
    }
};

/**
 * TODO description. maybe. already documented in geometryelement?
 * @private
 */
JXG.Line.prototype.update = function() {
    var i, funps;

    if(this.constrained) {
    	if(typeof this.funps != 'undefined') {
    		funps = this.funps();
    		this.point1 = funps[0];
    		this.point2 = funps[1];
    	} else {
            this.point1 = this.funp1();
            this.point2 = this.funp2();
    	}
    }

    if (this.needsUpdate) {
        if (true || !this.board.geonextCompatibilityMode) {
            this.updateStdform();
        }
        for(i=0; i<this.ticks.length; i++) {
            // i don't know why we need this, but if we don't check it, an error will be reported
            // when the origin is moved. it seems like this.ticks.length is lying.
            if(typeof this.ticks[i] != 'undefined')
                this.ticks[i].calculateTicksCoordinates();
        }
    }
    if(this.traced) {
        this.cloneToBackground(true);
    }
};

/**
 * TODO description. already documented in geometryelement?
 * @private
 */
JXG.Line.prototype.updateStdform = function() {
   /*
    var nx = -(this.point2.coords.usrCoords[2]-this.point1.coords.usrCoords[2]);
    var ny =  this.point2.coords.usrCoords[1]-this.point1.coords.usrCoords[1];
    var c = -(nx*this.point1.coords.usrCoords[1]+ny*this.point1.coords.usrCoords[2]);

    this.stdform[0] = c;
    this.stdform[1] = nx;
    this.stdform[2] = ny;
    */
    var v = JXG.Math.crossProduct(this.point1.coords.usrCoords,this.point2.coords.usrCoords);
    this.stdform[0] = v[0];
    this.stdform[1] = v[1];
    this.stdform[2] = v[2];
    this.stdform[3] = 0;
    this.normalize();
};

/**
 * Uses the boards renderer to update the line.
 * @private
 */
 JXG.Line.prototype.updateRenderer = function () {
    var wasReal;

    if (this.needsUpdate && this.visProp['visible']) {
        wasReal = this.isReal;
        this.isReal = (isNaN(this.point1.coords.usrCoords[1]+this.point1.coords.usrCoords[2]+this.point2.coords.usrCoords[1]+this.point2.coords.usrCoords[2]))?false:true;
        if (this.isReal) {
            if (wasReal!=this.isReal) {
                this.board.renderer.show(this);
                if(this.hasLabel && this.label.content.visProp['visible']) this.board.renderer.show(this.label.content);
            }
            this.board.renderer.updateLine(this);
        } else {
            if (wasReal!=this.isReal) {
                this.board.renderer.hide(this);
                if(this.hasLabel && this.label.content.visProp['visible']) this.board.renderer.hide(this.label.content);
            }
        }

        //this.board.renderer.updateLine(this); // Why should we need this?
        this.needsUpdate = false;
    }

    /* Update the label if visible. */
    if(this.hasLabel && this.label.content.visProp['visible'] && this.isReal) {
        //this.label.setCoordinates(this.coords);
        this.label.content.update();
        //this.board.renderer.updateLabel(this.label);
        this.board.renderer.updateText(this.label.content);
    }
};

/**
 * Used to generate a polynomial for a point p that lies on this line, i.e. p is collinear to {@link #point1}
 * and {@link #point2}.
 * @param p The point for that the polynomial is generated.
 * @return An array containing the generated polynomial.
 * @private
 */
JXG.Line.prototype.generatePolynomial = function (/** JXG.Point */ p) /** array */{
    var u1 = this.point1.symbolic.x,
        u2 = this.point1.symbolic.y,
        v1 = this.point2.symbolic.x,
        v2 = this.point2.symbolic.y,
        w1 = p.symbolic.x,
        w2 = p.symbolic.y;

    /*
     * The polynomial in this case is determined by three points being collinear:
     *
     *      U (u1,u2)      W (w1,w2)                V (v1,v2)
     *  ----x--------------x------------------------x----------------
     *
     *  The collinearity condition is
     *
     *      u2-w2       w2-v2
     *     -------  =  -------           (1)
     *      u1-w1       w1-v1
     *
     * Multiplying (1) with denominators and simplifying is
     *
     *    u2w1 - u2v1 + w2v1 - u1w2 + u1v2 - w1v2 = 0
     */

    return [['(',u2,')*(',w1,')-(',u2,')*(',v1,')+(',w2,')*(',v1,')-(',u1,')*(',w2,')+(',u1,')*(',v2,')-(',w1,')*(',v2,')'].join('')];
};

/**
 * Calculates the rise of the line.
 * @type float
 * @return The rise of the line.
 */
JXG.Line.prototype.getRise = function () {
    if (Math.abs(this.stdform[2])>=JXG.Math.eps) {
        return -this.stdform[0]/this.stdform[2];
    } else {
        return Infinity;
    }
};

/**
 * Calculates the slope of the line.
 * @type float
 * @return The slope of the line or Infinity if the line is parallel to the y-axis.
 */
JXG.Line.prototype.getSlope = function () {
    if (Math.abs(this.stdform[2])>=JXG.Math.eps) {
        return -this.stdform[1]/this.stdform[2];
    } else {
        return Infinity;
    }
};

/**
 * Determines whether the line is drawn beyond {@link #point1} and {@link #point2} and updates the line.
 * @param {boolean} straightFirst True if the Line shall be drawn beyond {@link #point1}, false otherwise.
 * @param {boolean} straightLast True if the Line shall be drawn beyond {@link #point2}, false otherwise.
 * @see #straightFirst
 * @see #straightLast
 * @private
 */
 JXG.Line.prototype.setStraight = function (straightFirst, straightLast) {
    this.visProp['straightFirst'] = straightFirst;
    this.visProp['straightLast'] = straightLast;

    this.board.renderer.updateLine(this);
};

/**
 * Determines whether the line has arrows at start or end of the line. Is stored in visProp['firstArrow'] and visProp['lastArrow']
 * @param {boolean} firstArrow True if there is an arrow at the start of the line, false otherwise.
 * @param {boolean} lastArrow True if there is an arrow at the end of the line, false otherwise.
 * @private
 */
JXG.Line.prototype.setArrow = function (firstArrow, lastArrow) {
     this.visProp['firstArrow'] = firstArrow;
     this.visProp['lastArrow'] = lastArrow;

     this.board.renderer.updateLine(this);
};

/**
 * Calculates TextAnchor. DESCRIPTION
 * @type JXG.Coords
 * @return Text anchor coordinates as JXG.Coords object.
 * @private
 */
JXG.Line.prototype.getTextAnchor = function() {
    return new JXG.Coords(JXG.COORDS_BY_USER, [0.5*(this.point2.X() - this.point1.X()),0.5*(this.point2.Y() - this.point1.Y())],this.board);
};

/**
 * Calculates LabelAnchor. DESCRIPTION
 * @type JXG.Coords
 * @return Text anchor coordinates as JXG.Coords object.
 * @private
 */
JXG.Line.prototype.getLabelAnchor = function() {
    var coords,screenCoords1,screenCoords2,
        relCoords, slope;

    if(!this.visProp['straightFirst'] && !this.visProp['straightLast']) {
        return new JXG.Coords(JXG.COORDS_BY_USER, [this.point2.X()-0.5*(this.point2.X() - this.point1.X()),this.point2.Y()-0.5*(this.point2.Y() - this.point1.Y())],this.board);
    }
    else {
        screenCoords1 = new JXG.Coords(JXG.COORDS_BY_USER, this.point1.coords.usrCoords, this.board);
        screenCoords2 = new JXG.Coords(JXG.COORDS_BY_USER, this.point2.coords.usrCoords, this.board);
        this.board.renderer.calcStraight(this, screenCoords1, screenCoords2);

        if(this.visProp['straightFirst']) {
            coords = screenCoords1;
        }
        else {
            coords = screenCoords2;
        }
        // Hack
        if(this.label.content != null) {
            relCoords = [0,0];
            slope = this.getSlope();
            if(coords.scrCoords[2]==0) {
                if(slope == Infinity) {
                    relCoords = [10,-10];
                }
                else if(slope >= 0) {
                    relCoords = [10,-10];
                }
                else {
                    relCoords = [-10,-10];
                }
            }
            else if(coords.scrCoords[2]==this.board.canvasHeight) {
                if(slope == Infinity) {
                    relCoords = [10,10];
                }
                else if(slope >= 0) {
                    relCoords = [-10,10];
                }
                else {
                    relCoords = [10,10];
                }
            }
            if(coords.scrCoords[1]==0) {
                if(slope == Infinity) {
                    relCoords = [10,10]; // ??
                }
                else if(slope >= 0) {
                    relCoords = [10,-10];
                }
                else {
                    relCoords = [10,10];
                }
            }
            else if(coords.scrCoords[1]==this.board.canvasWidth) {
                if(slope == Infinity) {
                    relCoords = [-10,10]; // ??
                }
                else if(slope >= 0) {
                    relCoords = [-10,10];
                }
                else {
                    relCoords = [-10,-10];
                }
            }
            this.label.content.relativeCoords = new JXG.Coords(JXG.COORDS_BY_USER, [relCoords[0]/this.board.stretchX,relCoords[1]/this.board.stretchY],this.board);
        }
        return coords;
    }
};

/**
 * Clone the element to the background to leave a trail of it on the board.
 * @param {boolean} addToTrace Not used.
 */
JXG.Line.prototype.cloneToBackground = function(addToTrace) {
    var copy = {},
        r, s;

    copy.id = this.id + 'T' + this.numTraces;
    this.numTraces++;
    copy.point1 = this.point1;
    copy.point2 = this.point2;

    copy.stdform = this.stdform;
    JXG.clearVisPropOld(copy);

    copy.board = {};
    copy.board.unitX = this.board.unitX;
    copy.board.unitY = this.board.unitY;
    copy.board.zoomX = this.board.zoomX;
    copy.board.zoomY = this.board.zoomY;
    copy.board.stretchX = this.board.stretchX;
    copy.board.stretchY = this.board.stretchY;
    copy.board.origin = this.board.origin;
    copy.board.canvasHeight = this.board.canvasHeight;
    copy.board.canvasWidth = this.board.canvasWidth;
    copy.board.dimension = this.board.dimension;
    copy.board.algebra = this.board.algebra;

    copy.visProp = this.visProp;
    s = this.getSlope();
    r = this.getRise();
    copy.getSlope = function() { return s; };
    copy.getRise = function() { return r; };

    this.board.renderer.enhancedRendering = true;
    this.board.renderer.drawLine(copy);
    this.board.renderer.enhancedRendering = false;
    this.traces[copy.id] = document.getElementById(copy.id);

    delete copy;

/*
    var id = this.id + 'T' + this.numTraces;
    this.traces[id] = this.board.renderer.cloneSubTree(this,id,'lines');
    this.numTraces++;
*/
};

/**
 * DESCRIPTION
 * @param transform A {@link #JXG.Transformation} object or an array of it.
 */
JXG.Line.prototype.addTransform = function (/** JXG.Transformation,array */ transform) {
    var list, i;
    if (JXG.isArray(transform)) {
        list = transform;
    } else {
        list = [transform];
    }
    for (i=0;i<list.length;i++) {
        this.point1.transformations.push(list[i]);
        this.point2.transformations.push(list[i]);
    }
};

/**
 * TODO DESCRIPTION. What is this method for? -- michael
 * @param method TYPE & DESCRIPTION. UNUSED.
 * @param x TYPE & DESCRIPTION
 * @param y TYPE & DESCRIPTION
 */
JXG.Line.prototype.setPosition = function (method, x, y) {
    //var oldCoords = this.coords;
    //if(this.group.length != 0) {
    // AW: Do we need this for lines?
        // this.coords = new JXG.Coords(method, [x,y], this.board);
        // this.group[this.group.length-1].dX = this.coords.scrCoords[1] - oldCoords.scrCoords[1];
        // this.group[this.group.length-1].dY = this.coords.scrCoords[2] - oldCoords.scrCoords[2];
        // this.group[this.group.length-1].update(this);
    //} else {
        var t = this.board.create('transform',[x,y],{type:'translate'});
        if (this.point1.transformations.length>0 && this.point1.transformations[this.point1.transformations.length-1].isNumericMatrix) {
            this.point1.transformations[this.point1.transformations.length-1].melt(t);
        } else {
            this.point1.addTransform(this.point1,t);
        }
        if (this.point2.transformations.length>0 && this.point2.transformations[this.point2.transformations.length-1].isNumericMatrix) {
            this.point2.transformations[this.point2.transformations.length-1].melt(t);
        } else {
            this.point2.addTransform(this.point2,t);
        }
        //this.addTransform(t);
        //this.update();
    //}
};

/**
 * Treat the line as parametric curve in homogeneuous coordinates.
 * <pre>x = 1 * sin(theta)*cos(phi)
 * y = 1 * sin(theta)*sin(phi)
 * z = 1 * sin(theta)</pre>
 * and the line is the set of solutions of <tt>a*x+b*y+c*z = 0</tt>.
 * It follows:
 * <pre>sin(theta)*(a*cos(phi)+b*sin(phi))+c*cos(theta) = 0</pre>
 * Define:
 * <pre>  A = (a*cos(phi)+b*sin(phi))
 *   B = c</pre>
 * Then
 * <pre>cos(theta) = A/sqrt(A*A+B*B)
 * sin(theta) = -B/sqrt(A*A+B*B)</pre>
 * and <tt>X(phi) = x</tt> from above.
 * phi runs from 0 to 1
 * @type float
 * @return X(phi) TODO description
 */
JXG.Line.prototype.X = function (phi) {
    var a = this.stdform[1],
        b = this.stdform[2],
        c = this.stdform[0],
        A, B, sq, sinTheta, cosTheta;
    phi *= Math.PI;
    A = a*Math.cos(phi)+b*Math.sin(phi);
    B = c;
    sq = Math.sqrt(A*A+B*B);
    sinTheta = -B/sq;
    cosTheta = A/sq;
    if (Math.abs(cosTheta)<this.board.algebra.eps) { cosTheta = 1.0; }
    return sinTheta*Math.cos(phi)/cosTheta;
};

/**
 * Treat the line as parametric curve in homogeneous coordinates. See {@link #X} for a detailed description.
 * @type float
 * @return Y(phi) TODO description
 */
JXG.Line.prototype.Y = function (phi) {
    var a = this.stdform[1],
        b = this.stdform[2],
        c = this.stdform[0],
        A, B, sq, sinTheta, cosTheta;
    phi *= Math.PI;
    A = a*Math.cos(phi)+b*Math.sin(phi);
    B = c;
    sq = Math.sqrt(A*A+B*B);
    sinTheta = -B/sq;
    cosTheta = A/sq;
    if (Math.abs(cosTheta)<this.board.algebra.eps) { cosTheta = 1.0; }
    return sinTheta*Math.sin(phi)/cosTheta;
};

/**
 * Treat the line as parametric curve in homogeneous coordinates. See {@link #X} for a detailed description.
 * @type float
 * @return Z(phi) TODO description
 */
JXG.Line.prototype.Z = function (phi) {
    var a = this.stdform[1],
        b = this.stdform[2],
        c = this.stdform[0],
        A, B, sq, cosTheta;
    phi *= Math.PI;
    A = a*Math.cos(phi)+b*Math.sin(phi);
    B = c;
    sq = Math.sqrt(A*A+B*B);
    cosTheta = A/sq;
    if (Math.abs(cosTheta)>=this.board.algebra.eps) {
        return 1.0;
    } else {
        return 0.0;
    }
};

/**
 * TODO circle?!? --michael
 * private or public? --michael
 * Treat the circle as parametric curve:
 * t runs from 0 to 1
 * @private
 */
JXG.Line.prototype.minX = function () {
    return 0.0;
};

/**
 * TODO circle?!? --michael
 * private or public? --michael
 * Treat the circle as parametric curve:
 * t runs from 0 to 1
 * @private
 */
JXG.Line.prototype.maxX = function () {
    return 1.0;
};

/**
 * Adds ticks to this line. Ticks can be added to any kind of line: line, arrow, and axis.
 * @param {JXG.Ticks} ticks Reference to a ticks object which is describing the ticks (color, distance, how many, etc.).
 * @type String
 * @return Id of the ticks object.
 */
JXG.Line.prototype.addTicks = function(ticks) {
    if(ticks.id == '' || typeof ticks.id == 'undefined')
        ticks.id = this.id + '_ticks_' + (this.ticks.length+1);

    this.board.renderer.drawTicks(ticks);
    this.ticks.push(ticks);

    this.ticks[this.ticks.length-1].updateRenderer();

    return ticks.id;
};

/**
 * Removes all ticks from a line.
 */
JXG.Line.prototype.removeAllTicks = function() {
    var t;
    for(t=this.ticks.length; t>0; t--) {
        this.board.renderer.remove(this.ticks[t-1].rendNode);
    }
    this.ticks = new Array();
};

/**
 * Removes ticks identified by parameter named tick from this line.
 * @param {JXG.Ticks} tick Reference to tick object to remove.
 */
JXG.Line.prototype.removeTicks = function(tick) {
    var t, j;
    if(this.defaultTicks != null && this.defaultTicks == tick) {
        this.defaultTicks = null;
    }

    for(t=this.ticks.length; t>0; t--) {
        if(this.ticks[t-1] == tick) {
            this.board.renderer.remove(this.ticks[t-1].rendNode);

            for(j=0; j<this.ticks[t-1].ticks.length; j++) {
                if(this.ticks[t-1].labels[j] != null)
                    if (this.ticks[t-1].labels[j].show) this.board.renderer.remove(this.ticks[t-1].labels[j].rendNode);
            }
            delete(this.ticks[t-1]);
        }
    }
};

/**
 * @class This element is used to provide a constructor for a general line. A general line is given by two points. By setting additional properties
 * a line can be used as an arrow and/or axis.
 * @pseudo
 * @description
 * @name Line
 * @augments JXG.Line
 * @constructor
 * @type JXG.Line
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point,array_JXG.Point,array} point1,point2 Parent elements can be two elements either of type {@link JXG.Point} or array of
 * numbers describing the coordinates of a point. In the latter case the point will be constructed automatically as a fixed invisible point.
 * @param {number_number_number} a,b,c A line can also be created providing three numbers. The line is then described by the set of solutions
 * of the equation <tt>a*x+b*y+c*z = 0</tt>.
 * @example
 * // Create a line using point and coordinates/
 * // The second point will be fixed and invisible.
 * var p1 = board.create('point', [4.5, 2.0]);
 * var l1 = board.create('line', [p1, [1.0, 1.0]]);
 * </pre><div id="c0ae3461-10c4-4d39-b9be-81d74759d122" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var glex1_board = JXG.JSXGraph.initBoard('c0ae3461-10c4-4d39-b9be-81d74759d122', {boundingbox: [-1, 7, 7, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var glex1_p1 = glex1_board.create('point', [4.5, 2.0]);
 *   var glex1_l1 = glex1_board.create('line', [glex1_p1, [1.0, 1.0]]);
 * </script><pre>
 * @example
 * // Create a point using three coordinates
 * var l1 = board.create('line', [1.0, -2.0, 3.0]);
 * </pre><div id="cf45e462-f964-4ba4-be3a-c9db94e2593f" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var glex2_board = JXG.JSXGraph.initBoard('cf45e462-f964-4ba4-be3a-c9db94e2593f', {boundingbox: [-1, 7, 7, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var glex2_l1 = glex2_board.create('line', [1.0, -2.0, 3.0]);
 * </script><pre>
 */
JXG.createLine = function(board, parents, atts) {
    var el, p1, p2, i,
        c = [];
        
    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'line','withLabel'), layer:null});

    var constrained = false;
    if (parents.length == 2) { // The line is defined by two points (or coordinates of two points)
        if (parents[0].length>1) { // point 1 given by coordinates
            p1 = board.create('point', parents[0], {visible:false,fixed:true});
        } else if (parents[0].elementClass == JXG.OBJECT_CLASS_POINT) {
            p1 =  JXG.getReference(board,parents[0]);
        } else if ((typeof parents[0] == 'function') && (parents[0]().elementClass == JXG.OBJECT_CLASS_POINT)) {
            p1 = parents[0]();
            constrained = true;
        } else
            throw new Error("JSXGraph: Can't create line with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");

        if (parents[1].length>1) { // point 2 given by coordinates
            p2 = board.create('point', parents[1], {visible:false,fixed:true});
        } else if (parents[1].elementClass == JXG.OBJECT_CLASS_POINT) {
            p2 =  JXG.getReference(board,parents[1]);
        } else if ((typeof parents[1] == 'function') && (parents[1]().elementClass == JXG.OBJECT_CLASS_POINT)) {
            p2 = parents[1]();
            constrained = true;
        } else
            throw new Error("JSXGraph: Can't create line with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
        el = new JXG.Line(board, p1.id, p2.id, atts['id'], atts['name'],atts['withLabel'],atts['layer']);
        if(constrained) {
        	el.constrained = true;
        	el.funp1 = parents[0];
        	el.funp2 = parents[1];
        }
    }
    else if (parents.length==3) {  // Line is defined by three coordinates
        // free line
        for (i=0;i<3;i++) {
            if (typeof parents[i]=='number') {
                c[i] = function(z){ return function() { return z; }; }(parents[i]);
            } else if (typeof parents[i]=='function') {
                c[i] = parents[i];
            } else {
                throw new Error("JSXGraph: Can't create line with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "' and '" + (typeof parents[2])+ "'.");
                return;
            }
        }
        // point 1: (0,c,-b)
        p1 = board.create('point',[
                function() { return 0.0;},
                function() { return c[2]();},
                function() { return -c[1]();}],{visible:false,name:' '});
        // point 2: (b^2+c^2,-ba+c,-ca-b)
        p2 = board.create('point',[
                function() { return c[2]()*c[2]()+c[1]()*c[1]();},
                function() { return -c[1]()*c[0]()+c[2]();},
                function() { return -c[2]()*c[0]()-c[1]();}],{visible:false,name:' '});
        el = new JXG.Line(board, p1.id, p2.id, atts['id'], atts['name'],atts['withLabel']);
    }
    else if ((parents.length==1) && (typeof parents[0] == 'function') && (parents[0]().length == 2) &&
    		 (parents[0]()[0].elementClass == JXG.OBJECT_CLASS_POINT) && (parents[0]()[1].elementClass == JXG.OBJECT_CLASS_POINT)) {
    	var ps = parents[0]();
        el = new JXG.Line(board, ps[0].id, ps[1].id, atts['id'], atts['name'],atts['withLabel'],atts['layer']);
        el.constrained = true;
        el.funps = parents[0];
    } else
        throw new Error("JSXGraph: Can't create line with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
    return el;
};

JXG.JSXGraph.registerElement('line', JXG.createLine);

/**
 * @class This element is used to provide a constructor for a segment. It's strictly spoken just a wrapper for element {@link Line} with {@link JXG.Line#straightFirst}
 * and {@link JXG.Line#straightLast} properties set to false.
 * @pseudo
 * @description
 * @name Segment
 * @augments JXG.Line
 * @constructor
 * @type JXG.Line
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point,array_JXG.Point,array} point1,point2 Parent elements can be two elements either of type {@link JXG.Point} or array of numbers describing the
 * coordinates of a point. In the latter case the point will be constructed automatically as a fixed invisible point.
 * @param {number_number_number} a,b,c A line can also be created providing three numbers. The line is then described by the set of solutions
 * of the equation <tt>a*x+b*y+c*z = 0</tt>.
 * @see Line
 * @example
 * // Create a segment providing two points.
 *   var p1 = board.create('point', [4.5, 2.0]);
 *   var p2 = board.create('point', [1.0, 1.0]);
 *   var l1 = board.create('segment', [p1, p2]);
 * </pre><div id="d70e6aac-7c93-4525-a94c-a1820fa38e2f" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var slex1_board = JXG.JSXGraph.initBoard('d70e6aac-7c93-4525-a94c-a1820fa38e2f', {boundingbox: [-1, 7, 7, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var slex1_p1 = slex1_board.create('point', [4.5, 2.0]);
 *   var slex1_p2 = slex1_board.create('point', [1.0, 1.0]);
 *   var slex1_l1 = slex1_board.create('segment', [slex1_p1, slex1_p2]);
 * </script><pre>
 */
 JXG.createSegment = function(board, parents, atts) {
    var el;

    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'line','withLabel'), layer:null});
    atts['straightFirst'] = false;
    atts['straightLast'] = false;
    el = board.create('line', parents, atts);

    return el;
};

JXG.JSXGraph.registerElement('segment', JXG.createSegment);

/**
 * @class This element is used to provide a constructor for arrow, which is just a wrapper for element {@link Line} with {@link JXG.Line#straightFirst}
 * and {@link JXG.Line#straightLast} properties set to false and {@link JXG.Line#lastArrow} set to true.
 * @pseudo
 * @description
 * @name Arrow
 * @augments JXG.Line
 * @constructor
 * @type JXG.Line
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point,array_JXG.Point,array} point1,point2 Parent elements can be two elements either of type {@link JXG.Point} or array of numbers describing the
 * coordinates of a point. In the latter case the point will be constructed automatically as a fixed invisible point.
 * @param {number_number_number} a,b,c A line can also be created providing three numbers. The line is then described by the set of solutions
 * of the equation <tt>a*x+b*y+c*z = 0</tt>.
 * @see Line
 * @example
 * // Create an arrow providing two points.
 *   var p1 = board.create('point', [4.5, 2.0]);
 *   var p2 = board.create('point', [1.0, 1.0]);
 *   var l1 = board.create('arrow', [p1, p2]);
 * </pre><div id="1d26bd22-7d6d-4018-b164-4c8bc8d22ccf" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var alex1_board = JXG.JSXGraph.initBoard('1d26bd22-7d6d-4018-b164-4c8bc8d22ccf', {boundingbox: [-1, 7, 7, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var alex1_p1 = alex1_board.create('point', [4.5, 2.0]);
 *   var alex1_p2 = alex1_board.create('point', [1.0, 1.0]);
 *   var alex1_l1 = alex1_board.create('arrow', [alex1_p1, alex1_p2]);
 * </script><pre>
 */
JXG.createArrow = function(board, parents, attributes) {
    var el;

    attributes = JXG.checkAttributes(attributes,{withLabel:JXG.readOption(board.options,'line','withLabel'), layer:null});
    //if ( (JXG.isPoint(parents[0])) && (JXG.isPoint(parents[1])) ) { // The constructability decision is delkegated to the line object
        el = board.create('line',parents,attributes);
        //el = new JXG.Line(board, parents[0], parents[1], attributes['id'], attributes['name'],attributes['withLabel']);
        el.setStraight(false,false);
        el.setArrow(false,true);
    //} // Ansonsten eine fette Exception um die Ohren hauen
    //else
    //    throw new Error("JSXGraph: Can't create arrow with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");

    return el;
};

JXG.JSXGraph.registerElement('arrow', JXG.createArrow);

/**
 * @class This element is used to provide a constructor for an axis. It's strictly spoken just a wrapper for element {@link Line} with {@link JXG.Line#straightFirst}
 * and {@link JXG.Line#straightLast} properties set to true. Additionally {@link JXG.Line#lastArrow} is set to true and default {@link Ticks} will be created.
 * @pseudo
 * @description
 * @name Axis
 * @augments JXG.Line
 * @constructor
 * @type JXG.Line
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point,array_JXG.Point,array} point1,point2 Parent elements can be two elements either of type {@link JXG.Point} or array of numbers describing the
 * coordinates of a point. In the latter case the point will be constructed automatically as a fixed invisible point.
 * @param {number_number_number} a,b,c A line can also be created providing three numbers. The line is then described by the set of solutions
 * of the equation <tt>a*x+b*y+c*z = 0</tt>.
 * @example
 * // Create an axis providing two coord pairs.
 *   var l1 = board.create('axis', [[0.0, 1.0], [1.0, 1.3]]);
 * </pre><div id="4f414733-624c-42e4-855c-11f5530383ae" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var axex1_board = JXG.JSXGraph.initBoard('4f414733-624c-42e4-855c-11f5530383ae', {boundingbox: [-1, 7, 7, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var axex1_l1 = axex1_board.create('axis', [[0.0, 1.0], [1.0, 1.3]]);
 * </script><pre>
 */
JXG.createAxis = function(board, parents, attributes) {
    var point1,
        point2,
        line, dist, c1, c2, len, defTicks;

    // Arrays oder Punkte, mehr brauchen wir nicht.
    if ( (JXG.isArray(parents[0]) || JXG.isPoint(parents[0]) ) && (JXG.isArray(parents[1]) || JXG.isPoint(parents[1])) ) {
        if( JXG.isPoint(parents[0]) )
            point1 = parents[0];
        else
            point1 = new JXG.Point(board, parents[0],'','',false);

        if( JXG.isPoint(parents[1]) )
            point2 = parents[1];
        else
            point2 = new JXG.Point(board,parents[1],'','',false);

        /* Make the points fixed */
        point1.fixed = true;
        point2.fixed = true;

        attributes = attributes || {};
        attributes.lastArrow = attributes.lastArrow || true;
        attributes.straightFirst = attributes.lastArrow || true;
        attributes.straightLast = attributes.straightLast || true;
        attributes.strokeWidth = attributes.strokeWidth || 1;
        attributes.withLabel = attributes.withLabel || false;
        attributes.highlightStrokeColor = attributes.highlightStrokeColor || attributes.strokeColor || board.options.axis.highlightStrokeColor;
        attributes.strokeColor = attributes.strokeColor || board.options.axis.strokeColor;

        line = board.create('line', [point1, point2], attributes);
        line.needsRegularUpdate = false;  // Axes only updated after zooming and moving of  the origin.

        attributes.minorTicks = attributes.minorTicks || 4;
        attributes.insertTicks = attributes.insertTicks || 'true';

        if(attributes.ticksDistance != 'undefined' && attributes.ticksDistance != null) {
            dist = attributes.ticksDistance;
        } else {
            c1 = new JXG.Coords(JXG.COORDS_BY_USER, [line.point1.coords.usrCoords.slice(1)],board);
            c2 = new JXG.Coords(JXG.COORDS_BY_USER, [line.point2.coords.usrCoords.slice(1)],board);
            board.renderer.calcStraight(line, c1, c2);
            len = c1.distance(JXG.COORDS_BY_USER,c2);
            //len *= 0.33;
            dist = 1.0; //len;
        }

        defTicks = board.create('ticks', [line, dist], attributes);
        defTicks.needsRegularUpdate = false;
        line.defaultTicks = defTicks;
    }
    else
        throw new Error("JSXGraph: Can't create point with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");

    return line;
};

JXG.JSXGraph.registerElement('axis', JXG.createAxis);

/**
 * @class With the element tangent the slope of a line, circle, or curve in a certain point can be visualized. A tangent is always constructed
 * by a glider on a line, circle, or curve and describes the tangent in the glider point on that line, circle, or curve.
 * @pseudo
 * @description
 * @name Tangent
 * @augments JXG.Line
 * @constructor
 * @type JXG.Line
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {Glider} g A glider on a line, circle, or curve.
 * @example
 * // Create a tangent providing a glider on a function graph
 *   var c1 = board.create('curve', [function(t){return t},function(t){return t*t*t;}]);
 *   var g1 = board.create('glider', [0.6, 1.2, c1]);
 *   var t1 = board.create('tangent', [g1]);
 * </pre><div id="7b7233a0-f363-47dd-9df5-4018d0d17a98" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var tlex1_board = JXG.JSXGraph.initBoard('7b7233a0-f363-47dd-9df5-4018d0d17a98', {boundingbox: [-6, 6, 6, -6], axis: true, showcopyright: false, shownavigation: false});
 *   var tlex1_c1 = tlex1_board.create('curve', [function(t){return t},function(t){return t*t*t;}]);
 *   var tlex1_g1 = tlex1_board.create('glider', [0.6, 1.2, tlex1_c1]);
 *   var tlex1_t1 = tlex1_board.create('tangent', [tlex1_g1]);
 * </script><pre>
 */
JXG.createTangent = function(board, parents, attributes) {
    var p,
        c,
        g, f, i, j, el, Dg, Df, tangent;

    if (parents.length==1) { // One arguments: glider on line, circle or curve
        p = parents[0];
        c = p.slideObject;
    } else if (parents.length==2) { // Two arguments: (point,line|curve|circle) or (line|curve|circle,point). // Not yet: curve!
        if (JXG.isPoint(parents[0])) {
            p = parents[0];
            c = parents[1];
        } else if (JXG.isPoint(parents[1])) {
            c = parents[0];
            p = parents[1];
        } else {
            throw new Error("JSXGraph: Can't create normal with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
        }
    } else {
        throw new Error("JSXGraph: Can't create normal with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
    }

    attributes = JXG.checkAttributes(attributes,{withLabel:JXG.readOption(board.options,'line','withLabel'), layer:null});
    
    if (c.elementClass == JXG.OBJECT_CLASS_LINE) {
        tangent = board.create('line', [c.point1,c.point2], attributes);
    } else if (c.elementClass == JXG.OBJECT_CLASS_CURVE) {
        if (c.curveType!='plot') {
            g = c.X;
            f = c.Y;
            tangent = board.create('line', [
                    function(){ return -p.X()*board.D(f)(p.position)+p.Y()*board.D(g)(p.position);},
                    function(){ return board.D(f)(p.position);},
                    function(){ return -board.D(g)(p.position);}
                    ], attributes );
            p.addChild(tangent);
            // this is required for the geogebra reader to display a slope
            tangent.glider = p;
        } else {  // curveType 'plot'
            // equation of the line segment: 0 = y*(x1-x2) + x*(y2-y1) + y1*x2-x1*y2
            tangent = board.create('line', [
                    function(){ i=Math.floor(p.position);
                                if (i==c.numberPoints-1) i--;
                                if (i<0) return 1.0;
                                return c.Y(i)*c.X(i+1)-c.X(i)*c.Y(i+1);},
                    function(){ i=Math.floor(p.position);
                                if (i==c.numberPoints-1) i--;
                                if (i<0) return 0.0;
                                return c.Y(i+1)-c.Y(i);},
                    function(){ i=Math.floor(p.position);
                                if (i==c.numberPoints-1) i--;
                                if (i<0) return 0.0;
                                return c.X(i)-c.X(i+1);}
                    ], attributes );
            p.addChild(tangent);
            // this is required for the geogebra reader to display a slope
            tangent.glider = p;
        }
    } else if (c.type == JXG.OBJECT_TYPE_TURTLE) {
            tangent = board.create('line', [
                    function(){ i=Math.floor(p.position);
                                for(j=0;j<c.objects.length;j++) {  // run through all curves of this turtle
                                    el = c.objects[j];
                                    if (el.type==JXG.OBJECT_TYPE_CURVE) {
                                        if (i<el.numberPoints) break;
                                        i-=el.numberPoints;
                                    }
                                }
                                if (i==el.numberPoints-1) i--;
                                if (i<0) return 1.0;
                                return el.Y(i)*el.X(i+1)-el.X(i)*el.Y(i+1);},
                    function(){ i=Math.floor(p.position);
                                for(j=0;j<c.objects.length;j++) {  // run through all curves of this turtle
                                    el = c.objects[j];
                                    if (el.type==JXG.OBJECT_TYPE_CURVE) {
                                        if (i<el.numberPoints) break;
                                        i-=el.numberPoints;
                                    }
                                }
                                if (i==el.numberPoints-1) i--;
                                if (i<0) return 0.0;
                                return el.Y(i+1)-el.Y(i);},
                    function(){ i=Math.floor(p.position);
                                for(j=0;j<c.objects.length;j++) {  // run through all curves of this turtle
                                    el = c.objects[j];
                                    if (el.type==JXG.OBJECT_TYPE_CURVE) {
                                        if (i<el.numberPoints) break;
                                        i-=el.numberPoints;
                                    }
                                }
                                if (i==el.numberPoints-1) i--;
                                if (i<0) return 0.0;
                                return el.X(i)-el.X(i+1);}
                    ], attributes );
            p.addChild(tangent);
            // this is required for the geogebra reader to display a slope
            tangent.glider = p;
    } else if (c.elementClass == JXG.OBJECT_CLASS_CIRCLE) {
        /*
        Dg = function(t){ return -c.Radius()*Math.sin(t); };
        Df = function(t){ return c.Radius()*Math.cos(t); };
        return board.create('line', [
                    function(){ return -p.X()*Df(p.position)+p.Y()*Dg(p.position);},
                    function(){ return Df(p.position);},
                    function(){ return -Dg(p.position);}
                    ], attributes );
        */
        // This construction should work on conics, too. p has to lie on c.
        board.create('line', [
                    function(){ return JXG.Math.matVecMult(c.quadraticform,p.coords.usrCoords)[0]; },
                    function(){ return JXG.Math.matVecMult(c.quadraticform,p.coords.usrCoords)[1]; },
                    function(){ return JXG.Math.matVecMult(c.quadraticform,p.coords.usrCoords)[2]; }
                ] , attributes);
        p.addChild(tangent);
        // this is required for the geogebra reader to display a slope
        tangent.glider = p;
    }
    
    return tangent;
};

/**
 * Register the element type tangent at JSXGraph
 * @private
 */
JXG.JSXGraph.registerElement('tangent', JXG.createTangent);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/
/** 
 * @fileoverview In this file the class Group is defined, a class for
 * managing grouping of points.
 * @author graphjs
 * @version 0.1
 */
 
/**
 * Creates a new instance of Group.
 * @class In this class all group management is done.
 * @param {String} id Unique identifier for this object.  If null or an empty string is given,
 * an unique id will be generated by Board
 * @param {String} name Not necessarily unique name, displayed on the board.  If null or an
 * empty string is given, an unique name will be generated.
 * @constructor
 */
JXG.Group = function(board, id, name) {
    var number,
        objArray,
        i, obj, el;
        
    this.board = board;
    this.objects = {};
    number = this.board.numObjects;
    this.board.numObjects++;

    if ((id == '') || (id == null) || (typeof id == 'undefined')) {
        this.id = this.board.id + 'Group' + number;
    } else {
        this.id = id;
    }
    
    this.type = JXG.OBJECT_TYPE_POINT;
    this.elementClass = JXG.OBJECT_CLASS_POINT;                

    if ((name == '') || (name == null) || (typeof name == 'undefined')) {
        this.name = 'group_' + this.board.generateName(this);
    } else {
        this.name = name;
    }
    delete(this.type);

    if ( (arguments.length == 4) && (JXG.isArray(arguments[3])) )
        objArray = arguments[3];
    else {
        objArray = [];
        for (i=3; i<arguments.length; i++) {
            objArray.push(arguments[i]);
        }
    }

    for (i=0; i<objArray.length; i++) {
        obj = JXG.getReference(this.board, objArray[i]);
        if( (!obj.fixed) && ( (obj.type == JXG.OBJECT_TYPE_POINT) || (obj.type == JXG.OBJECT_TYPE_GLIDER) ) ) {
            if (obj.group.length != 0) {
                this.addGroup(obj.group[obj.group.length-1]);
            } else {
                this.addPoint(obj);
            }
        }
    }
    
    for (el in this.objects) {
        this.objects[el].group.push(this);
    }

    this.dX = 0;
    this.dY = 0;
};

/**
 * Releases the group added to the points in this group, but only if this group is the last group.
 */
JXG.Group.prototype.ungroup = function() {
    var el;
    for (el in this.objects) {
        if (this.objects[el].group[this.objects[el].group.length-1] == this) {
            this.objects[el].group.pop();
        }
        delete(this.objects[el]);
    }
};

/**
 * Sends an update to all group members.
 * @param {JXG.Point} point The point that caused the update.
 */
JXG.Group.prototype.update = function(point) {
    var obj = null,
        el;
    
    for (el in this.objects) {
        obj = this.objects[el];
        if (obj.id != point.id) {
            obj.coords = new JXG.Coords(JXG.COORDS_BY_SCREEN, [obj.coords.scrCoords[1] + this.dX, obj.coords.scrCoords[2] + this.dY], obj.board);
        }
    }
    
    for (el in this.objects) {
        /* Wurde das Element vielleicht geloescht? */
        if (this.board.objects[el] != undefined) {
            /* Nein, wurde es nicht, also updaten */
            this.objects[el].update(false);
        } else { /* es wurde geloescht, also aus dem Array entfernen */
            delete(this.objects[el]);
        }
    }
};

/**
 * Adds an Point to this group.
 * @param {JXG.Point} object The object added to the group.
 */
JXG.Group.prototype.addPoint = function(object) {
    this.objects[object.id] = object;
};

/**
 * Adds an multiple points to this group.
 * @param {Array} objects An array of points to add to the group.
 */
JXG.Group.prototype.addPoints = function(objects) {
    var p;
    for (p in objects)
        this.objects[p.id] = p;
};

/**
 * Adds an Pint to this group.
 * @param {JXG.Point} object The object added to the group.
 */
JXG.Group.prototype.addGroup = function(group) {
    var el;
    for (el in group.objects) {
        this.addPoint(group.objects[el]);
    }
};

/**
 * Groups points.
 * @param {JXG.Board} board The board the points are on.
 * @param {Array} parents Array of points to group.
 * @param {Object} attributes Visual properties.
 * @type JXG.Group
 * @return An object of type JXG.Group.
 */
JXG.createGroup = function(board, parents, attributes) {
    return new JXG.Group(board, attributes["id"], attributes["name"], parents);
};

JXG.JSXGraph.registerElement('group', JXG.createGroup);
/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.

*/
    
/**
 * @fileoverview The geometry object Circle is defined in this file. Circle stores all
 * style and functional properties that are required to draw and move a circle on
 * a board.
 * @author graphjs
 * @version 0.1
 */

/**
 * A circle consists of all points with a given distance from one point. This point is called midpoint, the distance is called radius.
 * A circle can be constructed by providing a midpoint and a point on the circle or a midpoint and a radius (given as a number, function,
 * line, or circle). 
 * @class Creates a new circle object. Do not use this constructor to create a circle. Use {@link JXG.Board#create} with
 * type {@link Circle} instead.  
 * @constructor
 * @augments JXG.GeometryElement
 * @param {String,JXG.Board} board The board the new circle is drawn on.
 * @param {String} method Can be 
 * <ul><li> <b>'twoPoints'</b> which means the circle is defined by its midpoint and a point on the circle.</li>
 * <li><b>'pointRadius'</b> which means the circle is defined by its midpoint and its radius in user units</li>
 * <li><b>'pointLine'</b> which means the circle is defined by its midpoint and its radius given by the distance from the startpoint and the endpoint of the line</li>
 * <li><b>'pointCircle'</b> which means the circle is defined by its midpoint and its radius given by the radius of another circle</li></ul>
 * The parameters p1, p2 and radius must be set according to this method parameter.
 * @param {JXG.Point} p1 Midpoint of the circle.
 * @param {JXG.Point,JXG.Line,JXG.Circle} p2 Can be
 *<ul><li>a point on the circle if method is 'twoPoints'</li>
 <li>a line if the method is 'pointLine'</li>
 <li>a circle if the method is 'pointCircle'</li></ul>
 * @param {float} radius Only used when method is set to 'pointRadius'. Must be a given radius in user units.
 * @param {String} id Unique identifier for this object. If null or an empty string is given,
 * an unique id will be generated by Board
 * @param {String} name Not necessarily unique name. If null or an
 * empty string is given, an unique name will be generated.
 * @see JXG.Board#generateName
 */            

JXG.Circle = function (board, method, par1, par2, id, name, withLabel, layer) {
    /* Call the constructor of GeometryElement */
    this.constructor();

    /**
     * Type of element; Value is {@link JXG.OBJECT_TYPE_CIRCLE}.
     * @default {@link JXG.OBJECT_TYPE_CIRCLE}
     * @constant
     * @type number
     * @private
     */
    this.type = JXG.OBJECT_TYPE_CIRCLE;
    /**
     * Class of this element; Values is OBJECT_CLASS_CIRCLE.
     * @constant
     * @type number
     * @private
     */
    this.elementClass = JXG.OBJECT_CLASS_CIRCLE; 

    this.init(board, id, name);
    
    /**
     * Set the display layer.
     */
    if (layer == null) layer = board.options.layer['circle'];
    this.layer = layer;

    /**
     * Stores the given method.
     * Can be 
     * <ul><li><b>'twoPoints'</b> which means the circle is defined by its midpoint and a point on the circle.</li>
     * <li><b>'pointRadius'</b> which means the circle is defined by its midpoint and its radius given in user units or as term.</li>
     * <li><b>'pointLine'</b> which means the circle is defined by its midpoint and its radius given by the distance from the startpoint and the endpoint of the line.</li>
     * <li><b>'pointCircle'</b> which means the circle is defined by its midpoint and its radius given by the radius of another circle.</li></ul>
     * @type string
     * @see #midpoint
     * @see #point2
     * @see #radius
     * @see #line
     * @see #circle
     */
    this.method = method;
    
    /**
     * The circles midpoint. Do not set this parameter directly as it will break JSXGraph's update system.
     * @type JXG.Point
     */    
    this.midpoint = JXG.getReference(this.board, par1); 
    this.midpoint.addChild(this);
    
    /* documented in GeometryElement */
    this.visProp['visible'] = true;
    this.visProp['fillColor'] = this.board.options.circle.fillColor;
    this.visProp['highlightFillColor'] = this.board.options.circle.highlightFillColor;
    this.visProp['strokeColor'] = this.board.options.circle.strokeColor;
    this.visProp['highlightStrokeColor'] = this.board.options.circle.highlightStrokeColor;       
    
    /** Point on the circle only set if method equals 'twoPoints'. Do not set this parameter directly as it will break JSXGraph's update system.
     * @type JXG.Point
     * @see #method
     */
    this.point2 = null;
    
    /** Radius of the circle
     * only set if method equals 'pointRadius'
     * @type JXG.Point
     * @default null
     * @see #method     
     */    
    this.radius = 0;
    
    /** Line defining the radius of the circle given by the distance from the startpoint and the endpoint of the line
     * only set if method equals 'pointLine'. Do not set this parameter directly as it will break JSXGraph's update system.
     * @type JXG.Line
     * @default null
     * @see #method     
     */    
    this.line = null;
    
    /** Circle defining the radius of the circle given by the radius of the other circle
     * only set if method equals 'pointLine'. Do not set this parameter directly as it will break JSXGraph's update system.
     * @type JXG.Circle
     * @default null
     * @see #method     
     */     
    this.circle = null;

    if(method == 'twoPoints') {
        this.point2 = JXG.getReference(board,par2);
        this.point2.addChild(this);
        this.radius = this.Radius(); 
    }
    else if(method == 'pointRadius') {
        this.generateTerm(par2);  // Converts GEONExT syntax into JavaScript syntax
        this.updateRadius();                        // First evaluation of the graph  
    }
    else if(method == 'pointLine') {
        // dann ist p2 die Id eines Objekts vom Typ Line!
        this.line = JXG.getReference(board,par2);
        this.radius = this.line.point1.coords.distance(JXG.COORDS_BY_USER, this.line.point2.coords);    
    }
    else if(method == 'pointCircle') {
        // dann ist p2 die Id eines Objekts vom Typ Circle!
        this.circle = JXG.getReference(board,par2);
        this.radius = this.circle.Radius();     
    } 
    
    // create Label
    if (withLabel!=null) 
        this.createLabel(withLabel);
    
    if(method == 'twoPoints') {
        //this.point2 = JXG.getReference(board,par2);
        //this.point2.addChild(this);
        //this.radius = this.Radius(); 
        this.id = this.board.addCircle(this);           
    }
    else if(method == 'pointRadius') {
        //this.generateTerm(par2);  // Converts GEONExT syntax into JavaScript syntax
        //this.updateRadius();                        // First evaluation of the graph
        this.id = this.board.addCircle(this);
        this.notifyParents(par2);      
    }
    else if(method == 'pointLine') {
        // dann ist p2 die Id eines Objekts vom Typ Line!
        //this.line = JXG.getReference(board,par2);
        //this.radius = this.line.point1.coords.distance(JXG.COORDS_BY_USER, this.line.point2.coords);
        this.line.addChild(this);
        this.id = this.board.addCircle(this);        
    }
    else if(method == 'pointCircle') {
        // dann ist p2 die Id eines Objekts vom Typ Circle!
        //this.circle = JXG.getReference(board,par2);
        //this.radius = this.circle.Radius();
        this.circle.addChild(this);
        this.id = this.board.addCircle(this);        
    }    
};
JXG.Circle.prototype = new JXG.GeometryElement;

/**
 * Checks whether (x,y) is near the circle.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @return {bool} True if (x,y) is near the circle, False otherwise.
 * @private
 */
JXG.Circle.prototype.hasPoint = function (x, y) {
/*
    var genauigkeit = this.board.options.precision.hasPoint;
    genauigkeit = genauigkeit/(this.board.stretchX); 
    
    var checkPoint = new JXG.Coords(JXG.COORDS_BY_SCREEN, [x,y], this.board);
    var r = this.Radius();
    
    var dist = Math.sqrt(Math.pow(this.midpoint.coords.usrCoords[1]-checkPoint.usrCoords[1],2) + 
                         Math.pow(this.midpoint.coords.usrCoords[2]-checkPoint.usrCoords[2],2));
   
    return (Math.abs(dist-r) < genauigkeit);
*/    
    var prec = this.board.options.precision.hasPoint/(this.board.stretchX),
        mp = this.midpoint.coords.usrCoords,
        p = new JXG.Coords(JXG.COORDS_BY_SCREEN, [x,y], this.board),
        r = this.Radius();
    
    var dist = Math.sqrt((mp[1]-p.usrCoords[1])*(mp[1]-p.usrCoords[1]) + (mp[2]-p.usrCoords[2])*(mp[2]-p.usrCoords[2]));
    return (Math.abs(dist-r) < prec);
};

/**
 * Used to generate a polynomial for a point p that lies on this circle.
 * @param p The point for that the polynomial is generated.
 * @return An array containing the generated polynomial.
 * @private
 */
JXG.Circle.prototype.generatePolynomial = function (p) {
    /*
     * We have four methods to construct a circle:
     *   (a) Two points
     *   (b) Midpoint and radius
     *   (c) Midpoint and radius given by length of a segment
     *   (d) Midpoint and radius given by another circle
     *
     * In case (b) we have to distinguish two cases:
     *  (i)  radius is given as a number
     *  (ii) radius is given as a function
     * In the latter case there's no guarantee the radius depends on other geometry elements
     * in a polynomial way so this case has to be omitted.
     *
     * Another tricky case is case (d):
     * The radius depends on another circle so we have to cycle through the ancestors of each circle
     * until we reach one that's radius does not depend on another circles radius.
     *
     *
     * All cases (a) to (d) vary only in calculation of the radius. So the basic formulae for
     * a glider G (g1,g2) on a circle with midpoint M (m1,m2) and radius r is just:
     *
     *     (g1-m1)^2 + (g2-m2)^2 - r^2 = 0
     *
     * So the easiest case is (b) with a fixed radius given as a number. The other two cases (a)
     * and (c) are quite the same: Euclidean distance between two points A (a1,a2) and B (b1,b2),
     * squared:
     *
     *     r^2 = (a1-b1)^2 + (a2-b2)^2
     *
     * For case (d) we have to cycle recursively through all defining circles and finally return the
     * formulae for calculating r^2. For that we use JXG.Circle.symbolic.generateRadiusSquared().
     */

    var m1 = this.midpoint.symbolic.x;
    var m2 = this.midpoint.symbolic.y;
    var g1 = p.symbolic.x;
    var g2 = p.symbolic.y;

    var rsq = this.generateRadiusSquared();

    /* No radius can be calculated (Case b.ii) */
    if (rsq == '')
        return [];

    var poly = '((' + g1 + ')-(' + m1 + '))^2 + ((' + g2 + ')-(' + m2 + '))^2 - (' + rsq + ')';
    return [poly];
};

/**
 * Generate symbolic radius calculation for loci determination with Groebner-Basis algorithm.
 * @type String
 * @return String containing symbolic calculation of the circle's radius or an empty string
 * if the radius can't be expressed in a polynomial equation.
 * @private
 */
JXG.Circle.prototype.generateRadiusSquared = function () {
    /*
     * Four cases:
     *
     *   (a) Two points
     *   (b) Midpoint and radius
     *   (c) Midpoint and radius given by length of a segment
     *   (d) Midpoint and radius given by another circle
     */

    var rsq = '';

    if (this.method == "twoPoints") {
        var m1 = this.midpoint.symbolic.x;
        var m2 = this.midpoint.symbolic.y;
        var p1 = this.point2.symbolic.x;
        var p2 = this.point2.symbolic.y;

        rsq = '(' + p1 + '-' + m1 + ')^2 + (' + p2 + '-' + m2 + ')^2';
    } else if (this.method == "pointRadius") {
        if (typeof(this.radius) == 'number')
            rsq = '' + this.radius*this.radius;
    } else if (this.method == "pointLine") {
        var p1 = this.line.point1.symbolic.x;
        var p2 = this.line.point1.symbolic.y;

        var q1 = this.line.point2.symbolic.x;
        var q2 = this.line.point2.symbolic.y;

        rsq = '(' + p1 + '-' + q1 + ')^2 + (' + p2 + '-' + q2 + ')^2';
    } else if (this.method == "pointCircle") {
        rsq = this.circle.Radius();
    }

    return rsq;
};

/**
 * Uses the boards renderer to update the circle.
 */
JXG.Circle.prototype.update = function () {
    if(this.traced) {
        this.cloneToBackground(true);
    }
    
    if (this.needsUpdate) {
        if(this.method == 'pointLine') {
            this.radius = this.line.point1.coords.distance(JXG.COORDS_BY_USER, this.line.point2.coords); 
        }
        else if(this.method == 'pointCircle') {
            this.radius = this.circle.Radius();
        }
        else if(this.method == 'pointRadius') {
            this.radius = this.updateRadius();
        }
        if (!this.board.geonextCompatibilityMode) {
            this.updateStdform();
            this.updateQuadraticform();
        }
    }
};

/**
 * TODO description
 * @private
 */
JXG.Circle.prototype.updateQuadraticform = function () {
    var m = this.midpoint,
        mX = m.X(), mY = m.Y(), r = this.Radius();
    this.quadraticform = [[mX*mX+mY*mY-r*r,-mX,-mY],
                          [-mX,1,0],
                          [-mY,0,1]
                         ];
};

/**
 * TODO description
 * @private
 */
JXG.Circle.prototype.updateStdform = function () {
    this.stdform[3] = 0.5;
    this.stdform[4] = this.Radius();
    this.stdform[1] = -this.midpoint.coords.usrCoords[1];
    this.stdform[2] = -this.midpoint.coords.usrCoords[2];
    this.normalize();
};

/**
 * Uses the boards renderer to update the circle.
 * @private
 */
JXG.Circle.prototype.updateRenderer = function () {
/*
    if (this.needsUpdate) {
        this.board.renderer.updateCircle(this);
        this.needsUpdate = false;
    }
*/
    if (this.needsUpdate && this.visProp['visible']) {
        var wasReal = this.isReal;
        this.isReal = (isNaN(this.midpoint.coords.usrCoords[1]+this.midpoint.coords.usrCoords[2]+this.Radius()))?false:true;
        if (this.isReal) {
            if (wasReal!=this.isReal) { 
                this.board.renderer.show(this); 
                if(this.hasLabel && this.label.content.visProp['visible']) this.board.renderer.show(this.label.content); 
            }
            this.board.renderer.updateCircle(this);
        } else {
            if (wasReal!=this.isReal) { 
                this.board.renderer.hide(this); 
                if(this.hasLabel && this.label.content.visProp['visible']) this.board.renderer.hide(this.label.content); 
            }
        }
        this.needsUpdate = false;
    }
    
    /* Update the label if visible. */
    if(this.hasLabel && this.label.content.visProp['visible'] && this.isReal) {
        //this.label.setCoordinates(this.coords);
        this.label.content.update();
        //this.board.renderer.updateLabel(this.label);
        this.board.renderer.updateText(this.label.content);
    }    
};

/**
 * TODO description
 * @param term TODO type & description
 * @private
 */
JXG.Circle.prototype.generateTerm = function (term) {
    if (typeof term=='string') {
         var elements = this.board.elementsByName;
         // Convert GEONExT syntax into  JavaScript syntax
         var newTerm = this.board.algebra.geonext2JS(term+'');
         this.updateRadius = new Function('return ' + newTerm + ';');
    } else if (typeof term=='number') {
        this.updateRadius = function() { return term; };
    } else { // function
        this.updateRadius = term;
    }
};   

/**
 * TODO description
 * @param contentStr TODO type&description
 * @private
 */
JXG.Circle.prototype.notifyParents = function (contentStr) {
    var res = null;
    var elements = this.board.elementsByName;
    
    if (typeof contentStr == 'string') 
        this.board.algebra.findDependencies(this,contentStr+'');
};

/**
 * Calculates the radius of the circle.
 * @type float
 * @return The radius of the circle
 */
JXG.Circle.prototype.Radius = function() {
    if(this.method == 'twoPoints') {
        return(Math.sqrt(Math.pow(this.midpoint.coords.usrCoords[1]-this.point2.coords.usrCoords[1],2) + Math.pow(this.midpoint.coords.usrCoords[2]-this.point2.coords.usrCoords[2],2)));
    }
    else if(this.method == 'pointLine' || this.method == 'pointCircle') {
        return this.radius;
    }
    else if(this.method == 'pointRadius') {
        return this.updateRadius();
    }
};

/**
  * @deprecated
  */
JXG.Circle.prototype.getRadius = function() {
    return this.Radius();
};

/**
 * TODO description
 * @private
 */
JXG.Circle.prototype.getTextAnchor = function() {
    return this.midpoint.coords;
};

/**
 * TODO description
 * @private
 */
JXG.Circle.prototype.getLabelAnchor = function() {
    if(this.method == 'twoPoints') {
        var deltaX = this.midpoint.coords.usrCoords[1]-this.point2.coords.usrCoords[1];
        var deltaY = this.midpoint.coords.usrCoords[2]-this.point2.coords.usrCoords[2];
        return new JXG.Coords(JXG.COORDS_BY_USER, [this.midpoint.coords.usrCoords[1]+deltaX, this.midpoint.coords.usrCoords[2]+deltaY], this.board);
    }
    else if(this.method == 'pointLine' || this.method == 'pointCircle' || this.method == 'pointRadius') {
        return new JXG.Coords(JXG.COORDS_BY_USER, [this.midpoint.coords.usrCoords[1]-this.Radius(),this.midpoint.coords.usrCoords[2]], this.board);
    }
};


/**
 * Clone the circle to the background.
 * @param addToTrace Not used yet. Always true.
 */
JXG.Circle.prototype.cloneToBackground = function(/** boolean */ addToTrace) {
    var copy = {};
    copy.id = this.id + 'T' + this.numTraces;
    this.numTraces++;
    copy.midpoint = {};
    copy.midpoint.coords = this.midpoint.coords;
    var r = this.Radius();
    copy.Radius = function() { return r; };
    copy.getRadius = function() { return r; }; // deprecated
    
    copy.board = {};
    copy.board.unitX = this.board.unitX;
    copy.board.unitY = this.board.unitY;
    copy.board.zoomX = this.board.zoomX;
    copy.board.zoomY = this.board.zoomY;
    copy.board.stretchX = this.board.stretchX;
    copy.board.stretchY = this.board.stretchY;

    copy.visProp = this.visProp;
    JXG.clearVisPropOld(copy);
    
    this.board.renderer.drawCircle(copy);
    this.traces[copy.id] = document.getElementById(copy.id);

    delete copy;
};

/**
 * TODO description
 * @param transform TODO type&description
 * @private
 */
JXG.Circle.prototype.addTransform = function (transform) {
    var list;
    if (JXG.isArray(transform)) {
        list = transform;
    } else {
        list = [transform];
    }
    for (var i=0;i<list.length;i++) {
        this.midpoint.transformations.push(list[i]);
        if (this.method == 'twoPoints') {
            this.point2.transformations.push(list[i]);
        }
    }
};

/**
 * TODO description
 * @param method TODO
 * @param x TODO
 * @param y TODO
 * @private
 */
JXG.Circle.prototype.setPosition = function (method, x, y) {
    //if(this.group.length != 0) {
    // AW: Do we need this for lines?
    //} else {
    var t = this.board.create('transform',[x,y],{type:'translate'});
    this.addTransform(t);
        //this.update();
    //}
};

/**
* Treat the circle as parametric curve:
* Return <tt>X(t)= radius*cos(t)+centerX</tt>, where t runs from 0 to 1.
* @param t TODO description
* @return TODO description
*/
JXG.Circle.prototype.X = function (/** float */ t) /** float */ {
    t *= 2.0*Math.PI;
    return this.Radius()*Math.cos(t)+this.midpoint.coords.usrCoords[1];
};

/**
* Treat the circle as parametric curve:
* Return <tt>Y(t)= radius*cos(t)+centerX</tt>
* t runs from 0 to 1
* @param t TODO description
* @return TODO description
*/
JXG.Circle.prototype.Y = function (/** float */ t) /** float */ {
    t *= 2.0*Math.PI;
    return this.Radius()*Math.sin(t)+this.midpoint.coords.usrCoords[2];
};

/**
 * Treat the circle as parametric curve:
 * t runs from 0 to 1
 * TODO description
 * @private
 */
JXG.Circle.prototype.minX = function () {
    return 0.0;
};

/**
 * Treat the circle as parametric curve:
 * t runs from 0 to 1
 * TODO description
 * @private
 */
JXG.Circle.prototype.maxX = function () {
    return 1.0;
};

JXG.Circle.prototype.Area = function() {
    var r = this.Radius();
    return r*r*Math.PI;
};

/**
 * @class This element is used to provide a constructor for a circle. 
 * @pseudo
 * @description  A circle consists of all points with a given distance from one point. This point is called midpoint, the distance is called radius.
 * A circle can be constructed by providing a midpoint and a point on the circle or a midpoint and a radius (given as a number, function,
 * line, or circle). 
 * @name Circle
 * @augments JXG.Circle
 * @constructor
 * @type JXG.Circle
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point_number,JXG.Point,JXG.Line,JXG.Circle} midpoint,radius The midpoint must be given as a {@link JXG.Point}, but the radius can be given
 * as a number (which will create a circle with a fixed radius), another {@link JXG.Point}, a {@link JXG.Line} (the distance of start and end point of the
 * line will determine the radius), or another {@link JXG.Circle}.
 * @example
 * // Create a circle providing two points
 * var p1 = board.create('point', [2.0, 2.0]);
 * var p2 = board.create('point', [2.0, 0.0]);
 * var c1 = board.create('circle', [p1, p2]);
 * 
 * // Create another circle using the above circle
 * var p3 = board.create('point', [3.0, 2.0]);
 * var c2 = board.create('circle', [p3, c1]);
 * </pre><div id="5f304d31-ef20-4a8e-9c0e-ea1a2b6c79e0" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var cex1_board = JXG.JSXGraph.initBoard('5f304d31-ef20-4a8e-9c0e-ea1a2b6c79e0', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var cex1_p1 = cex1_board.create('point', [2.0, 2.0]);
 *   var cex1_p2 = cex1_board.create('point', [2.0, 0.0]);
 *   var cex1_c1 = cex1_board.create('circle', [cex1_p1, cex1_p2]);
 *   var cex1_p3 = cex1_board.create('point', [3.0, 2.0]);
 *   var cex1_c2 = cex1_board.create('circle', [cex1_p3, cex1_c1]);
 * </script><pre>
 */
JXG.createCircle = function(board, parentArr, atts) {
    var el, p, i;
    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'circle','withLabel'), layer:null});
    
    p = [];
    for (i=0;i<parentArr.length;i++) {
        if (JXG.isPoint(parentArr[i])) {
            p[i] = parentArr[i];              // Point
        } else if (parentArr[i].length>1) {
            p[i] = board.create('point', parentArr[i], {visible:false,fixed:true});  // Coordinates
        } else {
            p[i] = parentArr[i];              // Something else (number, function, string)
        }
    }
    if( parentArr.length==2 && JXG.isPoint(p[0]) && JXG.isPoint(p[1]) ) {
        // Point/Point
        el = new JXG.Circle(board, 'twoPoints', p[0], p[1], atts['id'], atts['name'],atts['withLabel'],atts['layer']);
    } else if( ( JXG.isNumber(p[0]) || JXG.isFunction(p[0]) || JXG.isString(p[0])) && JXG.isPoint(p[1]) ) {
        // Number/Point
        el = new JXG.Circle(board, 'pointRadius', p[1], p[0], atts['id'], atts['name'],atts['withLabel'],atts['layer']);
    } else if( ( JXG.isNumber(p[1]) || JXG.isFunction(p[1]) || JXG.isString(p[1])) && JXG.isPoint(p[0]) ) {
        // Point/Number
        el = new JXG.Circle(board, 'pointRadius', p[0], p[1], atts['id'], atts['name'],atts['withLabel'],atts['layer']);
    } else if( (p[0].type == JXG.OBJECT_TYPE_CIRCLE) && JXG.isPoint(p[1]) ) {
        // Circle/Point
        el = new JXG.Circle(board, 'pointCircle', p[1], p[0], atts['id'], atts['name'],atts['withLabel'],atts['layer']);
    } else if( (p[1].type == JXG.OBJECT_TYPE_CIRCLE) && JXG.isPoint(p[0])) {
        // Point/Circle
        el = new JXG.Circle(board, 'pointCircle', p[0], p[1], atts['id'], atts['name'],atts['withLabel'],atts['layer']);
    } else if( (p[0].type == JXG.OBJECT_TYPE_LINE) && JXG.isPoint(p[1])) {
        // Circle/Point
        el = new JXG.Circle(board, 'pointLine', p[1], p[0], atts['id'], atts['name'],atts['withLabel'],atts['layer']);
    } else if( (p[1].type == JXG.OBJECT_TYPE_LINE) && JXG.isPoint(p[0])) {
        // Point/Circle
        el = new JXG.Circle(board, 'pointLine', p[0], p[1], atts['id'], atts['name'],atts['withLabel'],atts['layer']);
    } else if( parentArr.length==3 && JXG.isPoint(p[0]) && JXG.isPoint(p[1]) && JXG.isPoint(p[2])) {
        // Circle through three points
        var arr = JXG.createCircumcircle(board, p, atts); // returns [center, circle]
        arr[0].setProperty({visible:false});
        return arr[1];
    } else
        throw new Error("JSXGraph: Can't create circle with parent types '" + (typeof parentArr[0]) + "' and '" + (typeof parentArr[1]) + "'.");
    
    return el;
};

JXG.JSXGraph.registerElement('circle', JXG.createCircle);

/*
    Copyright 2010
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph. If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview In this file the conic sections defined.
 */

/**
 * @class This element is used to provide a constructor for an ellipse. An ellipse is given by two points (the foci) and a third point on the the ellipse or 
 * the length of the major axis.
 * @pseudo
 * @description
 * @name Ellipse
 * @augments JXG.Curve
 * @constructor
 * @type JXG.Curve
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point,array_JXG.Point,array_JXG.Point,array} point1,point2,point3 Parent elements can be three elements either of type {@link JXG.Point} or array of
 * numbers describing the coordinates of a point. In the latter case the point will be constructed automatically as a fixed invisible point.
 * @param {JXG.Point,array_JXG.Point,array_number,function} point1,point2,number Parent elements can be two elements either of type {@link JXG.Point} or array of
 * numbers describing the coordinates of a point. The third parameter is a number/function which defines the length of the major axis
 * Optional parameters four and five are numbers which define the curve length (e.g. start/end). Default values are -pi and pi.
 * @example
 * // Create an Ellipse by three points
 * var A = board.create('point', [-1,4]);
 * var B = board.create('point', [-1,-4);
 * var C = board.create('point', [1,1]);
 * var el = board.create('ellipse',[A,B,C]);
 * </pre><div id="a4d7fb6f-8708-4e45-87f2-2379ae2bd2c0" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var glex1_board = JXG.JSXGraph.initBoard('a4d7fb6f-8708-4e45-87f2-2379ae2bd2c0', {boundingbox:[-6,6,6,-6], keepaspectratio:true, showcopyright: false, shownavigation: false});
 *   var A = glex1_board.create('point', [-1,4]);
 *   var B = glex1_board.create('point', [-1,-4);
 *   var C = glex1_board.create('point', [1,1]);
 *   var el = glex1_board.create('ellipse',[A,B,C]);
 * </script><pre>
 */
JXG.createEllipse = function(board, parents, atts) {
    var F = [],  // focus 1 and focus 2
        C, 
        majorAxis,
        i,
        rotationMatrix;

    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'conic','withLabel'), layer:null});
    // The foci and the third point are either points or coordinate arrays.
    for (i=0;i<2;i++) {
        if (parents[i].length>1) { // focus i given by coordinates
            F[i] = board.create('point', parents[i], {visible:false,fixed:true});
        } else if (JXG.isPoint(parents[i])) { // focus i given by point
            F[i] = JXG.getReference(board,parents[i]);
        } else if ((typeof parents[i] == 'function') && (parents[i]().elementClass == JXG.OBJECT_CLASS_POINT)) {  // given by function
            F[i] = parents[i]();
        } else if (JXG.isString(parents[i])) { // focus i given by point name
            F[i] = JXG.getReference(board,parents[i]);
        } else
            throw new Error("JSXGraph: Can't create Ellipse with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
    }
    if (JXG.isNumber(parents[2])) { // length of major axis
        majorAxis = JXG.createFunction(parents[2],board);
    } else {
        if (JXG.isPoint(parents[2])) {                                               // point on ellipse
            C = JXG.getReference(board,parents[2]);
        } else if (parents[2].length>1) {                                            // point on ellipse given by coordinates
            C = board.create('point', parents[2], {visible:false,fixed:true});
        } else if ((typeof parents[2] == 'function') && (parents[2]().elementClass == JXG.OBJECT_CLASS_POINT)) {  // given by function
            C = parents[2]();
        } else if (JXG.isString(parents[2])) {                                      // focus i given by point name
            C = JXG.getReference(board,parents[2]);
        } else {
            throw new Error("JSXGraph: Can't create Ellipse with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "' and '" + (typeof parents[2]) +"'.");
        }
        majorAxis = function(){ return C.Dist(F[0])+C.Dist(F[1]);};
    }

    if (typeof parents[4]=='undefined') parents[4] = 1.0001*Math.PI;   // to
    if (typeof parents[3]=='undefined') parents[3] = -1.0001*Math.PI;  // from

    atts = JXG.checkAttributes(atts,{curveType:'parameter'});

    var M = board.create('point', [
                function(){return (F[0].X()+F[1].X())*0.5;},
                function(){return (F[0].Y()+F[1].Y())*0.5;}
            ],{visible:false, name:'', withLabel:false});

    var transformFunc = function() {
            var ax = F[0].X(),
                ay = F[0].Y(),
                bx = F[1].X(),
                by = F[1].Y(),
                beta;

            // Rotate by the slope of the line [F[0],F[1]]
            var sgn = (bx-ax>0)?1:-1;
            if (Math.abs(bx-ax)>0.0000001) {
                beta = Math.atan2(by-ay,bx-ax)+ ((sgn<0)?Math.PI:0);
            } else {
                beta = ((by-ay>0)?0.5:-0.5)*Math.PI;
            }
            var m = [
                        [1,    0,             0],
                        [M.X(),Math.cos(beta),-Math.sin(beta)],
                        [M.Y(),Math.sin(beta), Math.cos(beta)]
                    ];
            return m;
        };

    var polarForm = function(phi,suspendUpdate) {
                var a = majorAxis()*0.5,
                    e = F[1].Dist(F[0])*0.5,
                    b = Math.sqrt(a*a-e*e);
                if (!suspendUpdate) {
                    rotationMatrix = transformFunc();
                }
                return JXG.Math.matVecMult(rotationMatrix,[1,a*Math.cos(phi),b*Math.sin(phi)]);
        };

    var curve = board.create('curve',
                    [function(phi,suspendUpdate) {return polarForm(phi,suspendUpdate)[1];},
                     function(phi,suspendUpdate) {return polarForm(phi,suspendUpdate)[2];},parents[3],parents[4]],atts);
    return curve;
};

/**
 * @class This element is used to provide a constructor for an hyperbola. An hyperbola is given by two points (the foci) and a third point on the the hyperbola or 
 * the length of the major axis.
 * @pseudo
 * @description
 * @name Hyperbola
 * @augments JXG.Curve
 * @constructor
 * @type JXG.Curve
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point,array_JXG.Point,array_JXG.Point,array} point1,point2,point3 Parent elements can be three elements either of type {@link JXG.Point} or array of
 * numbers describing the coordinates of a point. In the latter case the point will be constructed automatically as a fixed invisible point.
 * @param {JXG.Point,array_JXG.Point,array_number,function} point1,point2,number Parent elements can be two elements either of type {@link JXG.Point} or array of
 * numbers describing the coordinates of a point. The third parameter is a number/function which defines the length of the major axis
 * Optional parameters four and five are numbers which define the curve length (e.g. start/end). Default values are -pi and pi.
 * @example
 * // Create an Hyperbola by three points
 * var A = board.create('point', [-1,4]);
 * var B = board.create('point', [-1,-4);
 * var C = board.create('point', [1,1]);
 * var el = board.create('hyperbola',[A,B,C]);
 * </pre><div id="cf99049d-a3fe-407f-b936-27d76550f8c4" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var glex1_board = JXG.JSXGraph.initBoard('cf99049d-a3fe-407f-b936-27d76550f8c4', {boundingbox:[-6,6,6,-6], keepaspectratio:true, showcopyright: false, shownavigation: false});
 *   var A = glex1_board.create('point', [-1,4]);
 *   var B = glex1_board.create('point', [-1,-4);
 *   var C = glex1_board.create('point', [1,1]);
 *   var el = glex1_board.create('hyperbola',[A,B,C]);
 * </script><pre>
 */
JXG.createHyperbola = function(board, parents, atts) {
    var F = [],  // focus 1 and focus 2
        C, 
        majorAxis,
        i,
        rotationMatrix;

    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'conic','withLabel'), layer:null});
    // The foci and the third point are either points or coordinate arrays.
    for (i=0;i<2;i++) {
        if (parents[i].length>1) { // focus i given by coordinates
            F[i] = board.create('point', parents[i], {visible:false,fixed:true});
        } else if (JXG.isPoint(parents[i])) { // focus i given by point
            F[i] = JXG.getReference(board,parents[i]);
        } else if ((typeof parents[i] == 'function') && (parents[i]().elementClass == JXG.OBJECT_CLASS_POINT)) {  // given by function
            F[i] = parents[i]();
        } else if (JXG.isString(parents[i])) { // focus i given by point name
            F[i] = JXG.getReference(board,parents[i]);
        } else
            throw new Error("JSXGraph: Can't create Hyperbola with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
    }
    if (JXG.isNumber(parents[2])) { // length of major axis
        majorAxis = JXG.createFunction(parents[2],board);
    } else {
        if (JXG.isPoint(parents[2])) {                                               // point on ellipse
            C = JXG.getReference(board,parents[2]);
        } else if (parents[2].length>1) {                                            // point on ellipse given by coordinates
            C = board.create('point', parents[2], {visible:false,fixed:true});
        } else if ((typeof parents[2] == 'function') && (parents[2]().elementClass == JXG.OBJECT_CLASS_POINT)) {  // given by function
            C = parents[2]();
        } else if (JXG.isString(parents[2])) {                                      // focus i given by point name
            C = JXG.getReference(board,parents[2]);
        } else {
            throw new Error("JSXGraph: Can't create Hyperbola with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "' and '" + (typeof parents[2]) +"'.");
        }
        majorAxis = function(){ return C.Dist(F[0])-C.Dist(F[1]);};
    }

    if (typeof parents[4]=='undefined') parents[4] = 1.0001*Math.PI;   // to
    if (typeof parents[3]=='undefined') parents[3] = -1.0001*Math.PI;  // from

    atts = JXG.checkAttributes(atts,{curveType:'parameter'});

    var M = board.create('point', [
                function(){return (F[0].X()+F[1].X())*0.5;},
                function(){return (F[0].Y()+F[1].Y())*0.5;}
            ],{visible:false, name:'', withLabel:false});

    var transformFunc = function() {
            var ax = F[0].X(),
                ay = F[0].Y(),
                bx = F[1].X(),
                by = F[1].Y(),
                beta;

            // Rotate by the slope of the line [F[0],F[1]]
            var sgn = (bx-ax>0)?1:-1;
            if (Math.abs(bx-ax)>0.0000001) {
                beta = Math.atan2(by-ay,bx-ax)+ ((sgn<0)?Math.PI:0);
            } else {
                beta = ((by-ay>0)?0.5:-0.5)*Math.PI;
            }
            var m = [
                        [1,    0,             0],
                        [M.X(),Math.cos(beta),-Math.sin(beta)],
                        [M.Y(),Math.sin(beta), Math.cos(beta)]
                    ];
            return m;
        };

    /*
          * Hyperbola is defined by (a*sec(t),b*tan(t)) and sec(t) = 1/cos(t)
          */
    var polarForm = function(phi,suspendUpdate) {
                var a = majorAxis()*0.5,
                    e = F[1].Dist(F[0])*0.5,
                    b = Math.sqrt(-a*a+e*e);
                if (!suspendUpdate) {
                    rotationMatrix = transformFunc();
                }
                return JXG.Math.matVecMult(rotationMatrix,[1,a/Math.cos(phi),b*Math.tan(phi)]);
        };
    var curve = board.create('curve',
                    [function(phi,suspendUpdate) {return polarForm(phi,suspendUpdate)[1];},
                     function(phi,suspendUpdate) {return polarForm(phi,suspendUpdate)[2];},parents[3],parents[4]],atts);

    return curve;
};

/**
 * @class This element is used to provide a constructor for a parabola. A parabola is given by one point (the focus) and a line (the directrix).
 * @pseudo
 * @description
 * @name Parabola
 * @augments JXG.Curve
 * @constructor
 * @type JXG.Curve
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point,array_JXG.Line} point,line Parent elements are a point and a line.
 * Optional parameters three and four are numbers which define the curve length (e.g. start/end). Default values are -pi and pi.
 * @example
 * // Create a parabola by a point C and a line l.
 * var A = board.create('point', [-1,4]);
 * var B = board.create('point', [-1,-4);
 * var l = board.create('line', [A,B]);
 * var C = board.create('point', [1,1]);
 * var el = board.create('parabola',[C,l]);
 * </pre><div id="524d1aae-217d-44d4-ac58-a19c7ab1de36" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var glex1_board = JXG.JSXGraph.initBoard('524d1aae-217d-44d4-ac58-a19c7ab1de36', {boundingbox:[-6,6,6,-6], keepaspectratio:true, showcopyright: false, shownavigation: false});
 *   var A = glex1_board.create('point', [-1,4]);
 *   var B = glex1_board.create('point', [-1,-4);
 *   var l = glex1_board.create('line', [A,B]);
 *   var C = glex1_board.create('point', [1,1]);
 *   var el = glex1_board.create('parabola',[C,l]);
 * </script><pre>
 */
JXG.createParabola = function(board, parents, atts) {
    var F1 = parents[0], // focus
        l = parents[1],  // directrix
        rotationMatrix;

    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'conic','withLabel'), layer:null});
    if (parents[0].length>1) { // focus 1 given by coordinates
        F1 = board.create('point', parents[0], {visible:false,fixed:true});
    } else if (JXG.isPoint(parents[0])) { // focus i given by point
        F1 = JXG.getReference(board,parents[0]);
    } else if ((typeof parents[0] == 'function') && (parents[0]().elementClass == JXG.OBJECT_CLASS_POINT)) {  // given by function
        F1 = parents[0]();
    } else if (JXG.isString(parents[0])) { // focus i given by point name
        F1 = JXG.getReference(board,parents[0]);
    } else
        throw new Error("JSXGraph: Can't create Parabola with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");

    if (typeof parents[3]=='undefined') parents[3] = 10.0;   // to
    if (typeof parents[2]=='undefined') parents[2] = -10.0;  // from

    atts = JXG.checkAttributes(atts,{curveType:'parameter'});

    var M = board.create('point', [
                function() {
                    var v = [0,l.stdform[1],l.stdform[2]];
                    v = JXG.Math.crossProduct(v,F1.coords.usrCoords);
                    return board.algebra.meetLineLine(v,l.stdform,0).usrCoords;
                }
            ],{visible:false, name:'', withLabel:false});

    var transformFunc = function() {
            var beta = Math.atan(l.getSlope()),
                x = (M.X()+F1.X())*0.5,
                y = (M.Y()+F1.Y())*0.5;
            beta += (F1.Y()-M.Y()<0 || (F1.Y()==M.Y() && F1.X()>M.X()) ) ? Math.PI : 0;

            // Rotate by the slope of the line l (Leitlinie = directrix)
            var m = [
                        [1,    0,             0],
                        [x*(1-Math.cos(beta))+y*Math.sin(beta),Math.cos(beta),-Math.sin(beta)],
                        [y*(1-Math.cos(beta))-x*Math.sin(beta),Math.sin(beta), Math.cos(beta)]
                    ];
            return m;
        };

    var polarForm = function(t,suspendUpdate) {
                var e = M.Dist(F1)*0.5;
                if (!suspendUpdate) {
                    rotationMatrix = transformFunc();
                }
                return JXG.Math.matVecMult(rotationMatrix,[1,t+(M.X()+F1.X())*0.5,t*t/(e*4)+(M.Y()+F1.Y())*0.5]);
        };
    var curve = board.create('curve',
                    [function(t,suspendUpdate) {return polarForm(t,suspendUpdate)[1];},
                     function(t,suspendUpdate) {return polarForm(t,suspendUpdate)[2];},
                     parents[2],parents[3]],atts);

    return curve;
};

/**
 * 
 * @class This element is used to provide a constructor for a generic conic section uniquely defined by five points.
 * @pseudo
 * @description
 * @name Conic
 * @augments JXG.Curve
 * @constructor
 * @type JXG.Conic
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point,array_JXG.Point,array_JXG.Point,array_JXG.Point,array_JXG.Point,array_} point,point,point,point,point Parent elements are five points.
 * @param {number_number_number_number_number_number} 6 numbers (a_00,a_11,a_22,a_01,a_12,a_22)
 * @example
 * // Create a conic section through the points A, B, C, D, and E.
 *  var A = board.create('point', [1,5]);
 *  var B = board.create('point', [1,2]);
 *  var C = board.create('point', [2,0]);
 *  var D = board.create('point', [0,0]);
 *  var E = board.create('point', [-1,5]);
 *  var conic = board.create('conic',[A,B,C,D,E]);
 * </pre><div id="2d79bd6a-db9b-423c-9cba-2497f0b06320" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var glex1_board = JXG.JSXGraph.initBoard('2d79bd6a-db9b-423c-9cba-2497f0b06320', {boundingbox:[-6,6,6,-6], keepaspectratio:true, showcopyright: false, shownavigation: false});
 *   var A = glex1_board.create('point', [1,5]);
 *   var B = glex1_board.create('point', [1,2]);
 *   var C = glex1_board.create('point', [2,0]);
 *   var D = glex1_board.create('point', [0,0]);
 *   var E = glex1_board.create('point', [-1,5]);
 *   var conic = glex1_board.create('conic',[A,B,C,D,E]);
 * </script><pre>
 */
JXG.createConic = function(board, parents, atts) {
    var rotationMatrix, eigen, a, b, c, M,
        c1, c2, 
        points = [], i, definingMat, 
        givenByPoints = (parents.length==5)?true:false, 
        p = [];

    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'conic','withLabel'), layer:null});
    if (givenByPoints) {
        for (i=0;i<5;i++) {
            if (parents[i].length>1) { // point i given by coordinates
                points[i] = board.create('point', parents[i], {visible:false,fixed:true});
            } else if (JXG.isPoint(parents[i])) { // point i given by point
                points[i] = JXG.getReference(board,parents[i]);
            } else if ((typeof parents[i] == 'function') && (parents[i]().elementClass == JXG.OBJECT_CLASS_POINT)) {  // given by function
                points[i] = parents[i]();
            } else if (JXG.isString(parents[i])) { // point i given by point name
                points[i] = JXG.getReference(board,parents[i]);
            } else
                throw new Error("JSXGraph: Can't create Conic section with parent types '" + (typeof parents[i]) + "'.");
        }
    } else {
        /* Usual notation (x,y,z):
         *  [[A0,A3,A5],
         *   [A3,A1,A4],
         *   [A5,A4,A2]]. 
         * Our notation (z,x,y): 
         *  [[-A2   , A5*2.0, A4*0.5],
         *   [A5*2.0,    -A0, A3*0.5],
         *   [A4*0.5, A3*0.5,    -A1]] 
        */
        //definingMat = [[0,0,0],[0,0,0],[0,0,0]];
        //definingMat[0][0] = (JXG.isFunction(parent[0])) ? 
        M = [[-parents[2],parents[5]*2.0,parents[4]*0.5],
             [parents[5]*2.0,-parents[0],parents[3]*0.5],
             [parents[4]*0.5,parents[3]*0.5,-parents[1]]];
    }

    // sym(A) = A + A^t . Manipulates A in place.
    var sym = function(A) {
        var i, j;
        for (i=0;i<3;i++) {
            for (j=i;j<3;j++) {
                A[i][j] += A[j][i];
            }
        }
        for (i=0;i<3;i++) {
            for (j=0;j<i;j++) {
                A[i][j] = A[j][i];
            }
        }
        return A;
    };

    // degconic(v,w) = sym(v*w^t)
    var degconic = function(v,w) {
        var i, j, mat = [[0,0,0],[0,0,0],[0,0,0]];
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                mat[i][j] = v[i]*w[j];
            }
        }
        return sym(mat);
    };

    // (p^t*B*p)*A-(p^t*A*p)*B
    var fitConic = function(A,B,p)  {
        var pBp, pAp, Mv, M = [[0,0,0],[0,0,0],[0,0,0]], i, j;
        Mv = JXG.Math.matVecMult(B,p);
        pBp = JXG.Math.innerProduct(p,Mv);
        Mv = JXG.Math.matVecMult(A,p);
        pAp = JXG.Math.innerProduct(p,Mv);
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                M[i][j] = pBp*A[i][j]-pAp*B[i][j];
            }
        }
        return M;
    };

    var polarForm = function(phi,suspendUpdate) {
        var i, j, len, v;
        if (!suspendUpdate) {
            if (givenByPoints) {
                // Copy the point coordinate vectors
                for (i=0;i<5;i++) { 
                    p[i] = points[i].coords.usrCoords; 
                }
                // Compute the quadratic form
                c1 = degconic(JXG.Math.crossProduct(p[0],p[1]),JXG.Math.crossProduct(p[2],p[3]));
                c2 = degconic(JXG.Math.crossProduct(p[0],p[2]),JXG.Math.crossProduct(p[1],p[3]));
                M = fitConic(c1,c2,p[4]);
            }
            // Compute Eigenvalues and Eigenvectors
            eigen = JXG.Math.Numerics.Jacobi(M);
            // Scale the Eigenvalues such that the first Eigenvalue is positive
            if (eigen[0][0][0]<0) {
                eigen[0][0][0] *= (-1);
                eigen[0][1][1] *= (-1);
                eigen[0][2][2] *= (-1);
            }
            // Normalize the Eigenvectors
            for (i=0;i<3;i++) {
                len = 0.0;
                for (j=0;j<3;j++) {
                    len += eigen[1][j][i]*eigen[1][j][i];
                }
                len = Math.sqrt(len);
                for (j=0;j<3;j++) {
                    eigen[1][j][i] /= len;
                }
            }
            rotationMatrix = eigen[1];
            c = Math.sqrt(Math.abs(eigen[0][0][0]));
            a = Math.sqrt(Math.abs(eigen[0][1][1]));
            b = Math.sqrt(Math.abs(eigen[0][2][2]));
        }
        if (eigen[0][1][1]<0.0) {
            v = JXG.Math.matVecMult(rotationMatrix,[1/c,Math.cos(phi)/a,Math.sin(phi)/b]);
        } else if (eigen[0][2][2]<0.0) {
            v = JXG.Math.matVecMult(rotationMatrix,[Math.sin(phi)/c,Math.cos(phi)/a,1/b]);
        } 
        // Normalize
        v[1] /= v[0];
        v[2] /= v[0];
        v[0] = 1.0;
        return v;
    };

    var curve = board.create('curve',
            [function(phi,suspendUpdate) {return polarForm(phi,suspendUpdate)[1];},
             function(phi,suspendUpdate) {return polarForm(phi,suspendUpdate)[2];},
             -Math.PI,Math.PI],atts);
    return curve;
};

JXG.JSXGraph.registerElement('ellipse', JXG.createEllipse);
JXG.JSXGraph.registerElement('hyperbola', JXG.createHyperbola);
JXG.JSXGraph.registerElement('parabola', JXG.createParabola);
JXG.JSXGraph.registerElement('conic', JXG.createConic);


/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * Creates a new instance of Polygon.
 * @class Polygon stores all style and functional properties that are required
 * to draw a polygon on a board.
 * @param {JXG.Board} board Reference to the board the polygon is drawn on.
 * @param {Array} vertices Unique identifiers for the points defining the polygon.
 * Last point must be first point.
 * @param {Array} borders Unique identifiers for the derived borderlines of the polygon
 * @param {String} id Unique identifier for this object.  If null or an empty string is given,
 * an unique id will be generated by Board
 * @param {String} name Not necessarily unique name, displayed on the board.  If null or an
 * empty string is given, an unique name will be generated.
 * @see JXG.Board#addPolygon
 * @constructor
 * @extends JXG.GeometryElement
 */

JXG.Polygon = function (board, vertices, borders, id, name, withLines, withLabel, lineLabels, layer) {
    /* Call the constructor of GeometryElement */
    this.constructor();
    /**
     * Sets type of GeometryElement, value is OBJECT_TYPE_POLYGON.
     * @final
     * @type int
     */ 
    this.type = JXG.OBJECT_TYPE_POLYGON;
    this.elementClass = JXG.OBJECT_CLASS_AREA;                
    
    this.init(board, id, name);
    /**
     * Set the display layer.
     */
    if (layer == null) layer = board.options.layer['polygon'];
    this.layer = layer;
    
    if( (typeof withLines == 'undefined') || (withLines == null) ) {
        withLines = true;
    }
    if( (typeof lineLabels == 'undefined') || (lineLabels == null) ) {
        lineLabels = false;
    }    
        
    /**
     * Is the polygon bordered by lines?
     * @type bool
     */
    this.withLines = withLines;

    /**
     * References to the points defining the polygon.
     * Last vertex is the same as first vertex.
     * 
     * @type Array
     */    
    this.vertices = [];    
    for(var i=0; i<vertices.length; i++) {
       var vertex = JXG.getReference(this.board, vertices[i]);
       this.vertices[i] = vertex;
    }
    
    if((typeof borders == 'undefined') || (borders == null)) {
        borders = [];
        for(var i=0; i<vertices.length-1; i++) {
            borders[i] = {};
        }
    }
    
    if(this.vertices[this.vertices.length-1] != this.vertices[0]) {
        this.vertices.push(this.vertices[0]);
        borders.push({});
    }
    
    this.visProp['fillColor'] = this.board.options.polygon.fillColor;
    this.visProp['highlightFillColor'] = this.board.options.polygon.highlightFillColor;
    this.visProp['fillOpacity'] = this.board.options.polygon.fillOpacity;
    this.visProp['highlightFillOpacity'] = this.board.options.polygon.highlightFillOpacity;
    
    var l;
 
    /**
     * References to the borderlines of the polygon.
     * 
     * @type Array
     */  
    this.borders = [];
    if(withLines) {
        for(var i=0; i<this.vertices.length-1; i++) {
            /* create the borderlines */
            l = new JXG.Line(board, this.vertices[i], this.vertices[i+1], borders[i].id, borders[i].name, lineLabels, this.layer); // keine Labels?
            l.setStraight(false,false); // Strecke
            this.borders[i] = l;
            l.parentPolygon = this;
        }
    }
    
    /* Add polygon as child to defining points */
    for(var i=0; i<this.vertices.length-1; i++) { // last vertex is first vertex
        var vertex = JXG.getReference(this.board, this.vertices[i]);
        vertex.addChild(this);
    }
    
    // create label 
    this.createLabel(withLabel,[0,0]);
    
    /* Register polygon at board */
    this.id = this.board.addPolygon(this);
};
JXG.Polygon.prototype = new JXG.GeometryElement;

/**
 * Checks whether (x,y) is near the polygon.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @return {bool} Always false, because the polygons interior shall not be highlighted
 */
JXG.Polygon.prototype.hasPoint = function (x,y) {
    return false;
};

/**
 * Uses the boards renderer to update the polygon.
 */
JXG.Polygon.prototype.updateRenderer = function () {
    if (this.needsUpdate) {
        this.board.renderer.updatePolygon(this);
        this.needsUpdate = false;
    }
    if(this.hasLabel && this.label.content.visProp['visible']) {
        //this.label.setCoordinates(this.coords);
        this.label.content.update();
        //this.board.renderer.updateLabel(this.label);
        this.board.renderer.updateText(this.label.content);
    }    
};

/**
 * return TextAnchor
 */
JXG.Polygon.prototype.getTextAnchor = function() {
    var a = 0;
    var b = 0;
    var x = 0;
    var y = 0;
    a = x = this.vertices[0].X();
    b = y = this.vertices[0].Y();
    for (var i = 0; i < this.vertices.length; i++) {
        if (this.vertices[i].X() < a)
            a = this.vertices[i].X();
        if (this.vertices[i].X() > x)
            x = this.vertices[i].X();
        if (this.vertices[i].Y() > b)
            b = this.vertices[i].Y();
        if (this.vertices[i].Y() < y)
            y = this.vertices[i].Y();
    }
    return new JXG.Coords(JXG.COORDS_BY_USER, [(a + x)*0.5, (b + y)*0.5], this.board);
};

JXG.Polygon.prototype.getLabelAnchor = function() {
    var a = 0;
    var b = 0;
    var x = 0;
    var y = 0;
    a = x = this.vertices[0].X();
    b = y = this.vertices[0].Y();
    for (var i = 0; i < this.vertices.length; i++) {
        if (this.vertices[i].X() < a)
            a = this.vertices[i].X();
        if (this.vertices[i].X() > x)
            x = this.vertices[i].X();
        if (this.vertices[i].Y() > b)
            b = this.vertices[i].Y();
        if (this.vertices[i].Y() < y)
            y = this.vertices[i].Y();
    }
    return new JXG.Coords(JXG.COORDS_BY_USER, [(a + x)*0.5, (b + y)*0.5], this.board);
};

/**
 * Copy the element to the background.
 */
JXG.Polygon.prototype.cloneToBackground = function(addToTrace) {
    var copy = {};
    copy.id = this.id + 'T' + this.numTraces;
    this.numTraces++;
    copy.vertices = this.vertices;
    copy.visProp = this.visProp;
    JXG.clearVisPropOld(copy);
    
    this.board.renderer.drawPolygon(copy);

    this.traces[copy.id] = $(copy.id);

    delete copy;
};

JXG.createPolygon = function(board, parents, atts) {
    var el, i;

    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'polygon','withLabel'), layer:null});
    // Sind alles Punkte?
    for(i=0; i<parents.length; i++) {
        parents[i] = JXG.getReference(board, parents[i]);
        if(!JXG.isPoint(parents[i]))
            throw new Error("JSXGraph: Can't create polygon with parent types other than 'point'.");
    }
    
    el = new JXG.Polygon(board, parents, atts["borders"], atts["id"],atts["name"],atts["withLines"],
                        atts['withLabel'],atts['lineLabels'],atts['layer']);
    
    if(atts["withLines"] || true) {
    	for(i=0; i<el.borders.length; i++) {
    		el.borders[i].setProperty(atts);
    	}
    }

    return el;
};

JXG.JSXGraph.registerElement('polygon', JXG.createPolygon);

JXG.Polygon.prototype.hideElement = function() {
    this.visProp['visible'] = false;
    this.board.renderer.hide(this);

    if(this.withLines) {
        for(var i=0; i<this.borders.length; i++) {
            this.borders[i].hideElement();
        }
    }
    
    if (this.hasLabel && this.label!=null) {
        this.label.hiddenByParent = true;
        if(this.label.content.visProp['visible']) {
            this.board.renderer.hide(this.label.content);
        }
    }    
};

JXG.Polygon.prototype.showElement = function() {
    this.visProp['visible'] = true;
    this.board.renderer.show(this);

    if(this.withLines) {
        for(var i=0; i<this.borders.length; i++) {
            this.borders[i].showElement();
        }
    }
};

JXG.Polygon.prototype.Area = function() {
    //Surveyor's Formula
    var area=0, i;
    for(i=0; i<this.vertices.length-1; i++) {
        area += (this.vertices[i].X()*this.vertices[i+1].Y()-this.vertices[i+1].X()*this.vertices[i].Y()); // last vertex is first vertex
    }
    area /= 2.0;
    return Math.abs(area);
};

/**
 * @class Constructs a regular polygon. It needs two points which define the base line and the number of vertices.
 * @pseudo
 * @description Constructs a regular polygon. It needs two points which define the base line and the number of vertices, or a set of points.
 * @constructor
 * @name RegularPolygon
 * @type JXG.Polygon
 * @augments JXG.Polygon
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point_JXG.Point_Number} p1,p2,n The constructed regular polygon has n vertices and the base line defined by p1 and p2.
 * @example
 * var p1 = board.create('point', [0.0, 2.0]);
 * var p2 = board.create('point', [2.0, 1.0]);
 *
 * var pol = board.create('regularpolygon', [p1, p2, 5]);
 * </pre><div id="682069e9-9e2c-4f63-9b73-e26f8a2b2bb1" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var ccmex1_board = JXG.JSXGraph.initBoard('682069e9-9e2c-4f63-9b73-e26f8a2b2bb1', {boundingbox: [-1, 9, 9, -1], axis: false, showcopyright: false, shownavigation: false});
 *   var regpol_p1 = ccmex1_board.create('point', [0.0, 2.0]);
 *   var regpol_p2 = ccmex1_board.create('point', [2.0, 1.0]);
 *   var regpol_cc1 = ccmex1_board.create('regularpolygon', [ccmex1_p1, ccmex1_p2, 5]);
 * </script><pre>
 * @example
 * var p1 = board.create('point', [0.0, 2.0]);
 * var p2 = board.create('point', [0.0,-2.0]);
 * var p3 = board.create('point', [-2.0,0.0]);
 *
 * var pol = board.create('regularpolygon', [p1, p2, p3]);
 * </pre><div id="096a78b3-bd50-4bac-b958-3be5e7df17ed" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var ccmex2_board = JXG.JSXGraph.initBoard('096a78b3-bd50-4bac-b958-3be5e7df17ed', {boundingbox: [-1, 9, 9, -1], axis: false, showcopyright: false, shownavigation: false});
 *   var ccmex2_p1 = ccmex2_board.create('point', [0.0, 2.0]);
 *   var ccmex2_p2 = ccmex2_board.create('point', [0.0, -2.0]);
 *   var ccmex2_p3 = ccmex2_board.create('point', [-2.0,0.0]);
 *   var ccmex2_cc1 = ccmex2_board.create('regularpolygon', [ccmex2_p1, ccmex2_p2, ccmex2_p3]);
 * </script><pre>
 */
JXG.createRegularPolygon = function(board, parents, atts) {
    var el, i, n, p = [], rot, c, len, pointsExist;

    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'polygon','withLabel'), layer:null});
    if (JXG.isNumber(parents[parents.length-1]) && parents.length!=3) {
        throw new Error("JSXGraph: A regular polygon needs two point and a number as input.");
    }

    len = parents.length;
    n = parents[len-1];
    if ((!JXG.isNumber(n) && !JXG.isPoint(JXG.getReference(board, n))) || n<3) {
        throw new Error("JSXGraph: The third parameter has to be number greater than 2 or a point.");
    }
    
    if (JXG.isPoint(JXG.getReference(board, n))) {  // Regular polygon given by n points
        n = len;
        pointsExist = true;
    } else {
        len--;
        pointsExist = false;
    }
    // Sind alles Punkte? 
    for(i=0; i<len; i++) {
        parents[i] = JXG.getReference(board, parents[i]);
        if(!JXG.isPoint(parents[i]))
            throw new Error("JSXGraph: Can't create regular polygon if the first two parameters aren't points.");
    }

    p[0] = parents[0];
    p[1] = parents[1];
    for (i=2;i<n;i++) {
        rot = board.create('transform', [Math.PI*(2.0-(n-2)/n),p[i-1]], {type:'rotate'});
        if (pointsExist) {
            p[i] = parents[i];
            p[i].addTransform(parents[i-2],rot);
        } else {
            p[i] = board.create('point',[p[i-2],rot],{name:'', withLabel:false,fixed:true,face:'o',size:1});
        }
    }
    el = board.create('polygon',p,atts);

    return el;
};

JXG.JSXGraph.registerElement('regularpolygon', JXG.createRegularPolygon);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview In this file the geometry element Curve is defined.
 */

/**
 * Curves are the common object for function graphs, parametric curves, polar curves, adn data plots.
 * @class Creates a new curve object. Do not use this constructor to create a curve. Use {@link JXG.Board#create} with
 * type {@link Curve}, or {@link Functiongraph} instead.  
 * @augments JXG.GeometryElement
 * @param {string,JXG.Board} board The board the new curve is drawn on.
 * @param {Array} defining terms An array with the functon terms, data points of the curve.
 * @param {String} id Unique identifier for the point. If null or an empty string is given,
 *  an unique id will be generated by Board
 * @param {String} name Not necessarily unique name for the point. If null or an
 *  empty string is given, an unique name will be generated
 * @param {boolean} show False if the point is invisible, True otherwise
 * @see JXG.Board#generateName
 * @see JXG.Board#addCurve
  */
JXG.Curve = function (board, parents, id, name, withLabel, layer) {
    this.constructor();
 
    this.points = []; 

    this.type = JXG.OBJECT_TYPE_CURVE;
    this.elementClass = JXG.OBJECT_CLASS_CURVE;                
    
    this.init(board, id, name);

    /**
     * Set the display layer.
     */
    if (layer == null) layer = board.options.layer['curve'];
    this.layer = layer;

    /** Use the algorithm by Gillam and Hohenwarter for plotting.
      * If false the naive algorithm is used.
      * It is much slower, but the result is better.
      */
    this.doAdvancedPlot = this.board.options.curve.doAdvancedPlot;
    
    /** 
      * Number of points on curves after mouseUp, i.e. high quality output.
      * Only used if this.doAdvancedPlot==false
      * May be overwritten.
      **/
    this.numberPointsHigh = this.board.options.curve.numberPointsHigh;
    /** 
      * Number of points on curves after mousemove, i.e. low quality output.
      * Only used if this.doAdvancedPlot==false
      * May be overwritten.
      **/
    this.numberPointsLow = this.board.options.curve.numberPointsLow;
    /** 
      * Number of points on curves. This value changes
      * between numberPointsLow and numberPointsHigh.
      * It is set in {@link #updateCurve}.
      */
    this.numberPoints = this.numberPointsHigh; 

    this.visProp['strokeWidth'] = this.board.options.curve.strokeWidth;

    this.visProp['visible'] = true;
    this.dataX = null;
    this.dataY = null;

    /**
     * This is just for the hasPoint() method.
     * @type int
     */
    //this.r = this.board.options.precision.hasPoint;
    
    /**
     * The curveType is set in @see generateTerm and used in 
     * {@link updateCurve}
     * Possible values are:
     * 'none'
     * 'plot': Data plot
     * 'parameter': we can not distinguish function graphs and parameter curves
     * 'functiongraph': function graph
     * 'polar'
     * 'implicit' (not yet)
     *
     * Only parameter and plot are set directly.
     * polar is set with setProperties only.
     **/
    // this.curveType = 'none';
    this.curveType = null;

    if (parents[0]!=null) {
        this.varname = parents[0];
    } else {
        this.varname = 'x';
    }
    this.xterm = parents[1];  // function graphs: "x"
    this.yterm = parents[2];  // function graphs: e.g. "x^2"
    this.generateTerm(this.varname,this.xterm,this.yterm,parents[3],parents[4]);  // Converts GEONExT syntax into JavaScript syntax
    this.updateCurve();                        // First evaluation of the curve
    
    this.createLabel(withLabel);
    this.id = this.board.addCurve(this);
    
    if (typeof this.xterm=='string') {
        this.notifyParents(this.xterm);
    }
    if (typeof this.yterm=='string') {
        this.notifyParents(this.yterm);
    }
};
JXG.Curve.prototype = new JXG.GeometryElement;

/**
 * Gives the default value of the left bound for the curve.
 * May be overwritten in @see generateTerm.
 */
JXG.Curve.prototype.minX = function () {
    if (this.curveType=='polar') {
        return 0.0;
    } else {
        var leftCoords = new JXG.Coords(JXG.COORDS_BY_SCREEN, [0, 0], this.board);
        return leftCoords.usrCoords[1];
    }
};

/**
 * Gives the default value of the right bound for the curve.
 * May be overwritten in @see generateTerm.
 */
JXG.Curve.prototype.maxX = function () {
    var rightCoords;
    if (this.curveType=='polar') {
        return 2.0*Math.PI;
    } else {
        rightCoords = new JXG.Coords(JXG.COORDS_BY_SCREEN, [this.board.canvasWidth, 0], this.board);
        return rightCoords.usrCoords[1];
    }
};

/**
 * Checks whether (x,y) is near the curve.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @param {y} Find closest point on the curve to (x,y)
 * @return {bool} True if (x,y) is near the curve, False otherwise.
 */
JXG.Curve.prototype.hasPoint = function (x,y) {
    var t, dist = Infinity, 
        c, trans, i, j, tX, tY,
        xi, xi1, yi, yi1,
        lbda, x0, y0, x1, y1, xy, den,
        steps = this.numberPointsLow, 
        d = (this.maxX()-this.minX())/steps,
        prec = this.board.options.precision.hasPoint/(this.board.unitX*this.board.zoomX),
        checkPoint, len,
        suspendUpdate = true;

    prec = prec*prec;
    checkPoint = new JXG.Coords(JXG.COORDS_BY_SCREEN, [x,y], this.board);
    x = checkPoint.usrCoords[1];
    y = checkPoint.usrCoords[2];
    if (this.curveType=='parameter' || this.curveType=='polar' || this.curveType=='functiongraph') { 
        // Brute fore search for a point on the curve close to the mouse pointer
        len = this.transformations.length;
        for (i=0,t=this.minX(); i<steps; i++) {
            tX = this.X(t,suspendUpdate);
            tY = this.Y(t,suspendUpdate);
            for (j=0; j<len; j++) {
                trans = this.transformations[j];
                trans.update();
                c = JXG.Math.matVecMult(trans.matrix,[1,tX,tY]);
                tX = c[1];
                tY = c[2];
            }
            dist = (x-tX)*(x-tX)+(y-tY)*(y-tY);
            if (dist<prec) { return true; }
            t+=d;
        }  
    } else if (this.curveType == 'plot') {
        //$('debug').innerHTML +='. ';
        len = this.numberPoints; // Rough search quality
        for (i=0;i<len-1;i++) {
            xi = this.X(i);
            xi1 = this.X(i+1);
            //if (i!=xi) {
            //    yi = this.Y(xi);
            //    yi1 = this.Y(xi1);
            //} else {
                yi = this.Y(i);
                yi1 = this.Y(i+1);
               // $('debug').innerHTML = this.Y.toString();
            //}
            x1 = xi1 - xi;
            y1 = yi1-yi;
            
            x0 = x-xi; //this.X(i);
            y0 = y-yi; //this.Y(i);
            den = x1*x1+y1*y1;
            
            if (den>=JXG.Math.eps) {
                xy = x0*x1+y0*y1;
                lbda = xy/den;
                dist = x0*x0+y0*y0 - lbda*xy;
            } else {
                lbda = 0.0;
                dist = x0*x0+y0*y0;
            }
            if (lbda>=0.0 && lbda<=1.0 && dist<prec) { 
                return true; 
            } 
        }
        return false;
    } 
    return (dist<prec);
};

/**
  * Allocate points in the Coords array this.points
  */
JXG.Curve.prototype.allocatePoints = function () {
    var i, len;
    len = this.numberPoints;
    if (this.points.length<this.numberPoints) {
        for (i=this.points.length; i<len; i++) {
            this.points[i] = new JXG.Coords(JXG.COORDS_BY_USER, [0,0], this.board);
        }
    }
};

/**
 * Computes for equidistant points on the x-axis the values
 * of the function, {@link #updateCurve}
 * Then, the update function of the renderer
 * is called. 
 */
JXG.Curve.prototype.update = function () {
    if (this.needsUpdate) {
        this.updateCurve();
    }
    return this;
};

/**
 * Then, the update function of the renderer
 * is called. 
 */
JXG.Curve.prototype.updateRenderer = function () {
    if (this.needsUpdate) {
        this.board.renderer.updateCurve(this);
        this.needsUpdate = false;
    }
    
    /* Update the label if visible. */
    if(this.hasLabel && this.label.content.visProp['visible']) {
        //this.label.setCoordinates(this.coords);
        this.label.content.update();
        //this.board.renderer.updateLabel(this.label);
        this.board.renderer.updateText(this.label.content);
    }       
    return this;
};

/**
  * For dynamic dataplots updateCurve
  * can be used to compute new entries
  * for the arrays this.dataX and
  * this.dataY. It is used in @see updateCurve.
  * Default is an empty method, can be overwritten
  * by the user.
  */
JXG.Curve.prototype.updateDataArray = function () { return this; };

/**
 * Computes for equidistant points on the x-axis the values
 * of the function. @see #update
 * If the mousemove event triggers this update, we use only few
 * points. Otherwise, e.g. on mouseup, many points are used.
 */
JXG.Curve.prototype.updateCurve = function () {
    var len, mi, ma, x, y, i,
        suspendUpdate = false;
    
    this.updateDataArray();
    mi = this.minX();
    ma = this.maxX();

    // Discrete data points
    if (this.dataX!=null) { // x-coordinates are in an array
        this.numberPoints = this.dataX.length;
        len = this.numberPoints;
        this.allocatePoints();  // It is possible, that the array length has increased.
        for (i=0; i<len; i++) {
            x = i;
            if (this.dataY!=null) { // y-coordinates are in an array
                y = i;
            } else {
                y = this.X(x); // discrete x data, continuous y data
            }
            this.points[i].setCoordinates(JXG.COORDS_BY_USER, [this.X(x,suspendUpdate),this.Y(y,suspendUpdate)], false); // The last parameter prevents rounding in usr2screen().
            this.updateTransform(this.points[i]);
            suspendUpdate = true;
        }
    } else { // continuous x data
        if (this.doAdvancedPlot) {
            this.updateParametricCurve(mi,ma,len);
        } else {
            if (this.board.updateQuality==this.board.BOARD_QUALITY_HIGH) {
                this.numberPoints = this.numberPointsHigh;
            } else {
                this.numberPoints = this.numberPointsLow;
            }
            len = this.numberPoints;
            this.allocatePoints();  // It is possible, that the array length has increased.
            this.updateParametricCurveNaive(mi,ma,len);
        }
    }
    this.getLabelAnchor();
    return this;
};

JXG.Curve.prototype.updateParametricCurveNaive = function(mi,ma,len) {
    var i, t,
        suspendUpdate = false,
        stepSize = (ma-mi)/len;
        
    for (i=0; i<len; i++) {
        t = mi+i*stepSize;
        this.points[i].setCoordinates(JXG.COORDS_BY_USER, [this.X(t,suspendUpdate),this.Y(t,suspendUpdate)], false); // The last parameter prevents rounding in usr2screen().
        this.updateTransform(this.points[i]);
        suspendUpdate = true;
    }
    return this;
};

JXG.Curve.prototype.updateParametricCurve = function(mi,ma,len) {
    var i, t, t0,
        suspendUpdate = false,
        po = new JXG.Coords(JXG.COORDS_BY_USER, [0,0], this.board),
        x, y, x0, y0, top, depth,
        MAX_DEPTH,
        MAX_XDIST,
        MAX_YDIST,
        dyadicStack = [],
        depthStack = [],
        pointStack = [],
        divisors = [], 
        //xd_ = NaN, yd_ = NaN,
        distOK = false,
        j = 0;

    
    if (this.board.updateQuality==this.board.BOARD_QUALITY_LOW) {
        MAX_DEPTH = 12;
        MAX_XDIST = 12;
        MAX_YDIST = 12;
    } else {
        MAX_DEPTH = 20;
        MAX_XDIST = 2;
        MAX_YDIST = 2;
    }
    
    divisors[0] = ma-mi;
    for (i=1;i<MAX_DEPTH;i++) {
        divisors[i] = divisors[i-1]*0.5;
    }
    
    i = 1;
    dyadicStack[0] = 1;
    depthStack[0] = 0;
    t = mi;
    po.setCoordinates(JXG.COORDS_BY_USER, [this.X(t,suspendUpdate),this.Y(t,suspendUpdate)], false);
    suspendUpdate = true;
    x0 = po.scrCoords[1];
    y0 = po.scrCoords[2];
    t0 = t;
    
    t = ma;
    po.setCoordinates(JXG.COORDS_BY_USER, [this.X(t,suspendUpdate),this.Y(t,suspendUpdate)], false);
    x = po.scrCoords[1];
    y = po.scrCoords[2];
    
    pointStack[0] = [x,y];
    
    top = 1;
    depth = 0;

    this.points = [];
    this.points[j++] = new JXG.Coords(JXG.COORDS_BY_SCREEN, [x0, y0], this.board);
    
    do {
        distOK = this.isDistOK(x0,y0,x,y,MAX_XDIST,MAX_YDIST)||this.isSegmentOutside(x0,y0,x,y);
        while ( depth<MAX_DEPTH &&
               (!distOK || depth<3 /*|| (j>1 &&!this.bendOK(xd_,yd_,x-x0,y-y0))*/) &&
               !(!this.isSegmentDefined(x0,y0,x,y) && depth>8)
            ) {
            dyadicStack[top] = i;
            depthStack[top] = depth;
            pointStack[top] = [x,y];
            top++;
            
            i = 2*i-1;
            depth++;
            t = mi+i*divisors[depth];
            po.setCoordinates(JXG.COORDS_BY_USER, [this.X(t,suspendUpdate),this.Y(t,suspendUpdate)], false);
            x = po.scrCoords[1];
            y = po.scrCoords[2];
            distOK = this.isDistOK(x0,y0,x,y,MAX_XDIST,MAX_YDIST)||this.isSegmentOutside(x0,y0,x,y);
        }
        /*
        if (this.board.updateQuality==this.board.BOARD_QUALITY_HIGH && !this.isContinuous(t0,t,MAX_DEPTH)) {
            //$('debug').innerHTML += 'x ';
            this.points[j] = new JXG.Coords(JXG.COORDS_BY_SCREEN, [NaN, NaN], this.board);
            //this.points[j] = new JXG.Coords(JXG.COORDS_BY_SCREEN, [1, 1], this.board);
            j++;
        }
        */
        this.points[j] = new JXG.Coords(JXG.COORDS_BY_SCREEN, [x, y], this.board);
        this.updateTransform(this.points[j]);
        j++;
        //xd_ = x-x0;
        //yd_ = x-y0;
        x0 = x;
        y0 = y;
        t0 = t;
        
        top--;
        x = pointStack[top][0];
        y = pointStack[top][1];
        depth = depthStack[top]+1;
        i = dyadicStack[top]*2;
        
    } while (top != 0);
    this.numberPoints = this.points.length;
    //$('debug').innerHTML = ' '+this.numberPoints;
    return this;
        
};

JXG.Curve.prototype.isSegmentOutside = function (x0,y0,x1,y1) {
    if (y0<0 && y1<0) { return true; }
    else if (y0>this.board.canvasHeight && y1>this.board.canvasHeight) { return true; }
    else if (x0<0 && x1<0) { return true; }
    else if (x0>this.board.canvasWidth && x1>this.board.canvasWidth) { return true; }
    return false;
};

JXG.Curve.prototype.isDistOK = function (x0,y0,x1,y1,MAXX,MAXY) {
    if (isNaN(x0+y0+x1+y1)) { return false; }
    return (Math.abs(x1-x0)<MAXY && Math.abs(y1-y0)<MAXY);
};

JXG.Curve.prototype.isSegmentDefined = function (x0,y0,x1,y1) {
    if (isNaN(x0+y0) && isNaN(x1+y1)) { return false; }
    return true;
};
/*
JXG.Curve.prototype.isContinuous = function (t0, t1, MAX_ITER) {
    var left, middle, right, tm,
        iter = 0,
        initDist, dist = Infinity,
        dl, dr; 

    if (Math.abs(t0-t1)<JXG.Math.eps) { return true; }
    left = new JXG.Coords(JXG.COORDS_BY_USER, [0,0], this.board);
    middle = new JXG.Coords(JXG.COORDS_BY_USER, [0,0], this.board);
    right = new JXG.Coords(JXG.COORDS_BY_USER, [0,0], this.board);
    
    left.setCoordinates(JXG.COORDS_BY_USER, [this.X(t0,true),this.Y(t0,true)], false);
    right.setCoordinates(JXG.COORDS_BY_USER, [this.X(t1,true),this.Y(t1,true)], false);
    
    initDist = Math.max(Math.abs(left.scrCoords[1]-right.scrCoords[1]),Math.abs(left.scrCoords[2]-right.scrCoords[2]));
    while (iter++<MAX_ITER && dist>initDist*0.9) {
        tm = (t0+t1)*0.5;
        middle.setCoordinates(JXG.COORDS_BY_USER, [this.X(tm,true),this.Y(tm,true)], false);
        dl = Math.max(Math.abs(left.scrCoords[1]-middle.scrCoords[1]),Math.abs(left.scrCoords[2]-middle.scrCoords[2]));
        dr = Math.max(Math.abs(middle.scrCoords[1]-right.scrCoords[1]),Math.abs(middle.scrCoords[2]-right.scrCoords[2]));
        
        if (dl>dr) {
            dist = dl;
            t1 = tm;
        } else {
            dist = dr;
            t0 = tm;
        }
        if (Math.abs(t0-t1)<JXG.Math.eps) { return true;}
    }
    if (dist>initDist*0.9) {
        return false;
    } else {
        return true;
    }
};
*/

/*
JXG.Curve.prototype.bendOK = function (xd_,yd_,xd,yd) {
    var ip = xd_*xd+yd_*yd,
        MAX_BEND = Math.tan(45*Math.PI/180.0);

    if (isNaN(ip)) {
        return true;
    } else if (ip<=0.0) {
        return false;
    } else {
        return Math.abs(xd_*yd-yd_*xd)<MAX_BEND*ip;
    }
};
*/

JXG.Curve.prototype.updateTransform = function (p) {
    var t, c, i, 
        len = this.transformations.length;
    if (len==0) {
        return p;
    }
    for (i=0; i<len; i++) {
        t = this.transformations[i];
        t.update();
        c = JXG.Math.matVecMult(t.matrix,p.usrCoords);
        p.setCoordinates(JXG.COORDS_BY_USER,[c[1],c[2]]);
    }
    return p;
};

JXG.Curve.prototype.addTransform = function (transform) {
    var list, i, len;
    if (JXG.isArray(transform)) {
        list = transform;
    } else {
        list = [transform];
    }
    len = list.length;
    for (i=0; i<len; i++) {
        this.transformations.push(list[i]);
    }
    return this;
};

JXG.Curve.prototype.setPosition = function (method, x, y) {
    //if(this.group.length != 0) {
    // AW: Do we need this for lines?
    //} else {
    var t = this.board.create('transform',[x,y],{type:'translate'});
    if (this.transformations.length>0 && this.transformations[this.transformations.length-1].isNumericMatrix) {
        this.transformations[this.transformations.length-1].melt(t);
    } else {
        this.addTransform(t);
    }
    //this.update();
    //}
    return this;
};

/**
 * Converts the GEONExT syntax of the defining function term into JavaScript.
 * New methods X() and Y() for the Curve object are generated, further
 * new methods for minX() and maxX().
 *
 * Also, all objects whose name appears in the term are searched and
 * the curve is added as child to these objects. (Commented out!!!!)
 * @see Algebra
 * @see #geonext2JS.
 */
JXG.Curve.prototype.generateTerm = function (varname, xterm, yterm, mi, ma) {
    var fx, fy;

    // Generate the methods X() and Y()
    if (JXG.isArray(xterm)) {
        this.dataX = xterm;
        this.X = function(i) { return this.dataX[i]; };
        this.curveType = 'plot';
        this.numberPoints = this.dataX.length;
    } else {
        this.X = JXG.createFunction(xterm,this.board,varname);
        if (JXG.isString(xterm)) { 
            this.curveType = 'functiongraph'; 
        } else if (JXG.isFunction(xterm) || JXG.isNumber(xterm)) {
            this.curveType = 'parameter';
        }
    }

    if (JXG.isArray(yterm)) {
        this.dataY = yterm;
        this.Y = function(i) { 
                if (JXG.isFunction(this.dataY[i])) { 
                    return this.dataY[i](); 
                } else {
                    return this.dataY[i]; 
                }
            };
    } else {
        this.Y = JXG.createFunction(yterm,this.board,varname);
    }

    // polar form
    if (JXG.isFunction(xterm) && JXG.isArray(yterm)) {
        // Xoffset, Yoffset
        fx = JXG.createFunction(yterm[0],this.board,'');
        fy = JXG.createFunction(yterm[1],this.board,'');
        this.X = function(phi){return (xterm)(phi)*Math.cos(phi)+fx();};
        this.Y = function(phi){return (xterm)(phi)*Math.sin(phi)+fy();};
        this.curveType = 'polar';
    }

    // Set the bounds
    // lower bound
    if (mi!=null) this.minX = JXG.createFunction(mi,this.board,'');
    if (ma!=null) this.maxX = JXG.createFunction(ma,this.board,'');

/*    
    // Find dependencies
    var elements = this.board.elementsByName;
    for (el in elements) {
        if (el != this.name) {
            var s1 = "X(" + el + ")";
            var s2 = "Y(" + el + ")";
            if (xterm.indexOf(s1)>=0 || xterm.indexOf(s2)>=0 ||
                yterm.indexOf(s1)>=0 || yterm.indexOf(s2)>=0) {
                elements[el].addChild(this);
            }
        }
    }
*/    
};

/**
 * Finds dependencies in a given term and notifies the parents by adding the
 * dependent object to the found objects child elements.
 * @param {String} term String containing dependencies for the given object.
 */
JXG.Curve.prototype.notifyParents = function (contentStr) {
    //var res = null;
    //var elements = this.board.elementsByName;
    this.board.algebra.findDependencies(this,contentStr);
};

/**
 * Calculates LabelAnchor.
 * @type JXG.Coords
 * @return Text anchor coordinates as JXG.Coords object.
 */
JXG.Curve.prototype.getLabelAnchor = function() {
    var c = new JXG.Coords(JXG.COORDS_BY_SCREEN, [0, this.board.canvasHeight*0.5], this.board);
    c = this.board.algebra.projectCoordsToCurve(c.usrCoords[1],c.usrCoords[2],0.0,this)[0];
    return c;
};

/**
 * @class This element is used to provide a constructor for curve, which is just a wrapper for element {@link Curve}. 
 * A curve is a mapping from R to R^2. t mapsto (x(t),y(t)). The graph is drawn for t in the interval [a,b]. 
 * <p>
 * The following types of curves can be plotted:
 * <ul>
 *  <li> parametric curves: t mapsto (x(t),y(t)), where x() and y() are univariate functions.
 *  <li> polar curves: curves commonly written with polar equations like spirals and cardioids.
 *  <li> data plots: plot linbe segments through a given list of coordinates.
 * </ul>
 * @pseudo
 * @description
 * @name Curve
 * @augments JXG.Curve
 * @constructor
 * @type JXG.Curve
 *
 * @param {function,number_function,number_function,number_function,number} x,y,a_,b_ Parent elements for Parametric Curves. 
 *                     <p>
 *                     x describes the x-coordinate of the curve. It may be a function term in one variable, e.g. x(t). 
 *                     In case of x being of type number, x(t) is set to  a constant function.
 *                     this function at the values of the array.
 *                     <p>
 *                     y describes the y-coordinate of the curve. In case of a number, y(t) is set to the constant function
 *                     returning this number. 
 *                     <p>
 *                     Further parameters are an optional number or function for the left interval border a, 
 *                     and an optional number or function for the right interval border b. 
 *                     <p>
 *                     Default values are a=-10 and b=10.
 * @param {array_array,function,number} x,y Parent elements for Data Plots. 
 *                     <p>
 *                     x and y are arrays contining the x and y coordinates of the data points which are connected by
 *                     line segments. The individual entries of x and y may also be functions.
 *                     In case of x being an array the curve type is data plot, regardless of the second parameter and 
 *                     if additionally the second parameter y is a function term the data plot evaluates.
 * @param {function_array,function,number_function,number_function,number} r,offset_,a_,b_ Parent elements for Polar Curves. 
 *                     <p>
 *                     The first parameter is a function term r(phi) describing the polar curve.
 *                     <p>
 *                     The second parameter is the offset of the curve. It has to be
 *                     an array containing numbers or functions describing the offset. Default value is the origin [0,0].
 *                     <p>
 *                     Further parameters are an optional number or function for the left interval border a, 
 *                     and an optional number or function for the right interval border b. 
 *                     <p>
 *                     Default values are a=-10 and b=10.
 * @see JXG.Curve
 * @example
 * // Parametric curve
 * // Create a curve of the form (t-sin(t), 1-cos(t), i.e.
 * // the cycloid curve.
 *   var graph = board.create('curve', 
 *                        [function(t){ return t-Math.sin(t);}, 
 *                         function(t){ return 1-Math.cos(t);},
 *                         0, 2*Math.PI]
 *                     );
 * </pre><div id="af9f818b-f3b6-4c4d-8c4c-e4a4078b726d" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var c1_board = JXG.JSXGraph.initBoard('af9f818b-f3b6-4c4d-8c4c-e4a4078b726d', {boundingbox: [-1, 5, 7, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var graph1 = c1_board.create('curve', [function(t){ return t-Math.sin(t);},function(t){ return 1-Math.cos(t);},0, 2*Math.PI]);
 * </script><pre>
 * @example
 * // Data plots
 * // Connect a set of points given by coordinates with dashed line segments.
 * // The x- and y-coordinates of the points are given in two separate 
 * // arrays.
 *   var x = [0,1,2,3,4,5,6,7,8,9];
 *   var y = [9.2,1.3,7.2,-1.2,4.0,5.3,0.2,6.5,1.1,0.0];
 *   var graph = board.create('curve', [x,y], {dash:2});
 * </pre><div id="7dcbb00e-b6ff-481d-b4a8-887f5d8c6a83" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var c3_board = JXG.JSXGraph.initBoard('7dcbb00e-b6ff-481d-b4a8-887f5d8c6a83', {boundingbox: [-1,10,10,-1], axis: true, showcopyright: false, shownavigation: false});
 *   var x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
 *   var y = [9.2, 1.3, 7.2, -1.2, 4.0, 5.3, 0.2, 6.5, 1.1, 0.0];
 *   var graph3 = c3_board.create('curve', [x,y], {dash:2});
 * </script><pre>
 * @example
 * // Polar plot
 * // Create a curve with the equation r(phi)= a*(1+phi), i.e.
 * // a cardioid.
 *   var a = board.create('slider',[[0,2],[2,2],[0,1,2]]);
 *   var graph = board.create('curve', 
 *                        [function(phi){ return a.Value()*(1-Math.cos(phi));}, 
 *                         [1,0], 
 *                         0, 2*Math.PI]
 *                     );
 * </pre><div id="d0bc7a2a-8124-45ca-a6e7-142321a8f8c2" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var c2_board = JXG.JSXGraph.initBoard('d0bc7a2a-8124-45ca-a6e7-142321a8f8c2', {boundingbox: [-3,3,3,-3], axis: true, showcopyright: false, shownavigation: false});
 *   var a = c2_board.create('slider',[[0,2],[2,2],[0,1,2]]);
 *   var graph2 = c2_board.create('curve', [function(phi){ return a.Value()*(1-Math.cos(phi));}, [1,0], 0, 2*Math.PI]);
 * </script><pre>
 */
JXG.createCurve = function(board, parents, attributes) {
    attributes = JXG.checkAttributes(attributes,{withLabel:JXG.readOption(board.options,'curve','withLabel'), layer:null});
    return new JXG.Curve(board, ['x'].concat(parents), attributes['id'], attributes['name'],   
                         attributes['withLabel'],attributes['layer']);
};

JXG.JSXGraph.registerElement('curve', JXG.createCurve);

/**
 * @class This element is used to provide a constructor for functiongraph, which is just a wrapper for element {@link Curve} with {@link JXG.Curve#X()}
 * set to x. The graph is drawn for x in the interval [a,b].
 * @pseudo
 * @description
 * @name Functiongraph
 * @augments JXG.Curve
 * @constructor
 * @type JXG.Curve
 * @param {function_number,function_number,function} f,a_,b_ Parent elements are a function term f(x) describing the function graph. 
 *         <p>
 *         Further, an optional number or function for the left interval border a, 
 *         and an optional number or function for the right interval border b. 
 *         <p>
 *         Default values are a=-10 and b=10.
 * @see JXG.Curve
 * @example
 * // Create a function graph for f(x) = 0.5*x*x-2*x
 *   var graph = board.create('functiongraph', 
 *                        [function(x){ return 0.5*x*x-2*x;}, -2, 4]
 *                     );
 * </pre><div id="efd432b5-23a3-4846-ac5b-b471e668b437" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var alex1_board = JXG.JSXGraph.initBoard('efd432b5-23a3-4846-ac5b-b471e668b437', {boundingbox: [-3, 7, 5, -3], axis: true, showcopyright: false, shownavigation: false});
 *   var graph = alex1_board.create('functiongraph', [function(x){ return 0.5*x*x-2*x;}, -2, 4]);
 * </script><pre>
 * @example
 * // Create a function graph for f(x) = 0.5*x*x-2*x with variable interval
 *   var s = board.create('slider',[[0,4],[3,4],[-2,4,5]]);
 *   var graph = board.create('functiongraph', 
 *                        [function(x){ return 0.5*x*x-2*x;}, 
 *                         -2, 
 *                         function(){return s.Value();}]
 *                     );
 * </pre><div id="4a203a84-bde5-4371-ad56-44619690bb50" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var alex2_board = JXG.JSXGraph.initBoard('4a203a84-bde5-4371-ad56-44619690bb50', {boundingbox: [-3, 7, 5, -3], axis: true, showcopyright: false, shownavigation: false});
 *   var s = alex2_board.create('slider',[[0,4],[3,4],[-2,4,5]]);
 *   var graph = alex2_board.create('functiongraph', [function(x){ return 0.5*x*x-2*x;}, -2, function(){return s.Value();}]);
 * </script><pre>
 */
JXG.createFunctiongraph = function(board, parents, attributes) {
    var par = ["x","x"].concat(parents);
    attributes = JXG.checkAttributes(attributes,{withLabel:JXG.readOption(board.options,'curve','withLabel'), layer:null});
    attributes['curveType'] = 'functiongraph';
    return new JXG.Curve(board, par, attributes['id'], attributes['name'],attributes['withLabel'],attributes['layer']);
};

JXG.JSXGraph.registerElement('functiongraph', JXG.createFunctiongraph);


/**
 * TODO
 * Create a dynamic spline interpolated curve given by sample points p_1 to p_n.
 * @param {JXG.Board} board Reference to the board the spline is drawn on.
 * @param {Array} parents Array of points the spline interpolates
 * @param {Object} attributes Define color, width, ... of the spline
 * @type JXG.Curve
 * @return Returns reference to an object of type JXG.Curve.
 */
JXG.createSpline = function(board, parents, attributes) {
    var F;
    attributes = JXG.checkAttributes(attributes,{withLabel:JXG.readOption(board.options,'curve','withLabel'), layer:null});
    F = function() {
        var D, x=[], y=[];
        
        var fct = function (t,suspended) {
            var i, j;
        
            if (!suspended) {
                x = [];
                y = [];

                // given as [x[], y[]]
                if(parents.length == 2 && JXG.isArray(parents[0]) && JXG.isArray(parents[1]) && parents[0].length == parents[1].length) {
                    for(i=0; i<parents[0].length; i++) {
                        if(typeof parents[0][i] == 'function')
                            x.push(parents[0][i]());
                        else
                            x.push(parents[0][i]);
                        if(typeof parents[1][i] == 'function')
                            y.push(parents[1][i]());
                        else
                            y.push(parents[1][i]);
                    }
                } else {
                    for(i=0; i<parents.length; i++) {
                        if(JXG.isPoint(parents[i])) {
                            //throw new Error("JSXGraph: JXG.createSpline: Parents has to be an array of JXG.Point.");
                            x.push(parents[i].X());
                            y.push(parents[i].Y());
                        } else if (JXG.isArray(parents[i]) && parents[i].length == 2) {     // given as [[x1,y1], [x2, y2], ...]
                            for(i=0; i<parents.length; i++) {
                                if(typeof parents[i][0] == 'function')
                                    x.push(parents[i][0]());
                                else
                                    x.push(parents[i][0]);
                                if(typeof parents[i][1] == 'function')
                                    y.push(parents[i][1]());
                                else
                                    y.push(parents[i][1]);
                            }
                        }
                    }
                }
        
                // The array D has only to be calculated when the position of one or more sample point
                // changes. otherwise D is always the same for all points on the spline.
                D = JXG.Math.Numerics.splineDef(x, y);
            }
            return JXG.Math.Numerics.splineEval(t, x, y, D);
        };
        return fct;
    };
    return new JXG.Curve(board, ["x","x", F()], attributes["id"], attributes["name"],
                        attributes['withLabel'],attributes['layer']);
};

/**
 * Register the element type spline at JSXGraph
 * @private
 */
JXG.JSXGraph.registerElement('spline', JXG.createSpline);

/**
 * @class This element is used to provide a constructor for Riemann sums, which is relaized as a special curve. 
 * @pseudo
 * @description
 * @name Riemannsum
 * @augments JXG.Curve
 * @constructor
 * @type JXG.Curve
 * @param {function_number,function_string,function_function,number_function,number} f,n,type_,a_,b_ Parent elements of Riemannsum are a 
 *         function term f(x) describing the function graph which is filled by the Riemann rectangles.
 *         <p>
 *         n determines the number of rectangles, it is either a fixed number or a function.
 *         <p>
 *         type is a string or function returning one of the values:  'left', 'right', 'middle', 'lower', 'upper', or 'trapezodial'.
 *         Default value is 'left'.
 *         <p>
 *         Further parameters are an optional number or function for the left interval border a, 
 *         and an optional number or function for the right interval border b. 
 *         <p>
 *         Default values are a=-10 and b=10.
 * @see JXG.Curve
 * @example
 * // Create Riemann sums for f(x) = 0.5*x*x-2*x.
 *   var s = board.create('slider',[[0,4],[3,4],[0,4,10]],{snapWidth:1});
 *   var f = function(x) { return 0.5*x*x-2*x; };
 *   var r = board.create('riemannsum', 
 *               [f, function(){return s.Value();}, 'upper', -2, 5],
 *               {fillOpacity:0.4}
 *               );
 *   var g = board.create('functiongraph',[f, -2, 5]);
 * </pre><div id="940f40cc-2015-420d-9191-c5d83de988cf" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var rs1_board = JXG.JSXGraph.initBoard('940f40cc-2015-420d-9191-c5d83de988cf', {boundingbox: [-3, 7, 5, -3], axis: true, showcopyright: false, shownavigation: false});
 *   var f = function(x) { return 0.5*x*x-2*x; };
 *   var s = rs1_board.create('slider',[[0,4],[3,4],[0,4,10]],{snapWidth:1});
 *   var r = rs1_board.create('riemannsum', [f, function(){return s.Value();}, 'upper', -2, 5], {fillOpacity:0.4});
 *   var g = rs1_board.create('functiongraph', [f, -2, 5]);
 * </script><pre>
 */
JXG.createRiemannsum = function(board, parents, attributes) {
    var n, type, f, par, c;
    
    attributes = JXG.checkAttributes(attributes,
                    {withLabel:JXG.readOption(board.options,'curve','withLabel'),layer:null,fillOpacity:0.3,fillColor:'#ffff00', curveType:'plot'});

    f = parents[0]; 
    n = JXG.createFunction(parents[1],board,'');
    if (n==null) {
        throw new Error("JSXGraph: JXG.createRiemannsum: argument '2' n has to be number or function.");
    }
    type = JXG.createFunction(parents[2],board,'',false);
    if (type==null) {
        throw new Error("JSXGraph: JXG.createRiemannsum: argument 3 'type' has to be string or function.");
    }

    par = ['x', [0], [0]].concat(parents.slice(3));
    /**
     * @private
     */
    c = new JXG.Curve(board, par, attributes['id'], attributes['name'], attributes['withLabel'],attributes['layer']);
    /**
     * @private
     */
    c.updateDataArray = function() {
            var u = JXG.Math.Numerics.riemann(f,n(),type(),this.minX(),this.maxX());
            this.dataX = u[0];
            this.dataY = u[1];
        };
    return c;
};

JXG.JSXGraph.registerElement('riemannsum', JXG.createRiemannsum);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/** 
 * @fileoverview In this file the geometry object Arc is defined. Arc stores all
 * style and functional properties that are required to draw an arc on a board.
 * @author graphjs
 * @version 0.1
 */
 
/**
 * Creates a new instance of Arc.
 * @class Arc stores all style and functional properties that are required
 * to draw an arc on a board.
 * @param {JXG.Board} board Reference to the board the arc is drawn on.
 * @param {JXG.Point} p1 Midpoint of the arc.
 * @param {JXG.Point} p2 Point defining the arcs radius
 * @param {JXG.Point} p3 This point defines the angle of the arcs section.
 * @param {String} id Unique identifier for this object.  If null or an empty string is given,
 * an unique id will be generated by Board
 * @param {String} name Not necessarily unique name, displayed on the board.  If null or an
 * empty string is given, an unique name will be generated.
 * @see JXG.Board#addArc
 * @constructor
 * @extends JXG.GeometryElement
 */
JXG.Arc = function (board, p1, p2, p3, id, name, withLabel, layer) {
    /* Call the constructor of GeometryElement */
    this.constructor();
    
    /**
     * Type of GeometryElement, value is OBJECT_TYPE_ARC.
     * @final
     * @type int
     */
    this.type = JXG.OBJECT_TYPE_ARC;
    
    /**
     * Class of the element, value is OBJECT_CLASS_CIRCLE.
     * @final
     * @type int
     */
    this.elementClass = JXG.OBJECT_CLASS_CIRCLE;

    /* Call init defined in GeometryElement to set board, id and name property */
    this.init(board, id, name);
    
    /**
     * Set the display layer.
     */
    if (layer == null) layer = board.options.layer['arc'];
    this.layer = layer;

    /**
     * Midpoint of the arc.
     * @type JXG.Point
     */
    this.midpoint = JXG.getReference(this.board, p1);
    /**
     * Point defining the arcs circle.
     * @type JXG.Point
     */
    this.point2 = JXG.getReference(this.board, p2);
    /**
     * The point defining the angle of the arc.
     * @type JXG.Point
     */
    this.point3 = JXG.getReference(this.board, p3);
    
    /**
     * This is just for the hasPoint() method. Precision for highlighting.
     * @type int
     */
    //this.r = this.board.options.precision.hasPoint;

    this.visProp['visible'] = true;
    
    this.visProp['firstArrow'] = this.board.options.arc.firstArrow;
    this.visProp['lastArrow'] = this.board.options.arc.lastArrow;
    
    this.visProp['fillColor'] = this.board.options.arc.fillColor;
    this.visProp['highlightFillColor'] = this.board.options.arc.highlightFillColor;
    this.visProp['strokeColor'] = this.board.options.arc.strokeColor;
    this.visProp['highlightStrokeColor'] = this.board.options.arc.highlightStrokeColor;     

    // create Label
    this.createLabel(withLabel,[0,0]);
    
    /* Register arc at board. */
    this.id = this.board.addArc(this);
    
    /* Add arc as child to defining points */
    this.midpoint.addChild(this);
    this.point2.addChild(this);
    this.point3.addChild(this);
};

JXG.Arc.prototype = new JXG.GeometryElement;

/**
 * Checks whether (x,y) is near the arc.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @return {bool} True if (x,y) is near the arc, False otherwise.
 */
JXG.Arc.prototype.hasPoint = function (x, y) { 
    var genauigkeit = this.board.options.precision.hasPoint/(this.board.stretchX);
    
    var checkPoint = new JXG.Coords(JXG.COORDS_BY_SCREEN, [x,y], this.board);
    var r = this.Radius();
    
    var dist = Math.sqrt(Math.pow(this.midpoint.coords.usrCoords[1]-checkPoint.usrCoords[1],2) + 
                         Math.pow(this.midpoint.coords.usrCoords[2]-checkPoint.usrCoords[2],2));
   
    var has = (Math.abs(dist-r) < genauigkeit);
    if(has) {
        var p = {};
        p.coords = new JXG.Coords(JXG.COORDS_BY_USER, 
                              [this.midpoint.coords.usrCoords[1], 
                               this.board.origin.usrCoords[2]/(this.board.stretchY)],
                              this.board);
        var angle1 = this.board.algebra.trueAngle(this.point2, this.midpoint, p);
        var angle2 = this.board.algebra.trueAngle(this.point3, this.midpoint, p);

        var xy = {};
        xy.coords = checkPoint;
        var angle3 = this.board.algebra.trueAngle(xy, this.midpoint, p); 
        if(angle1 >= angle2) {
            if(angle1 < angle3 || angle3 < angle2) {
                has = false;
            }
        }
        else {
            if(angle3 > angle1) {
                if(angle3 < angle2) {
                    has = false;
                }
            }
        }
    }
    return has;    
};

/**
 * Checks whether (x,y) is within the sector defined by the arc.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @return {bool} True if (x,y) is within the sector defined by the arc, False otherwise.
 */
JXG.Arc.prototype.hasPointSector = function (x, y) { 
    var genauigkeit = this.board.options.precision.hasPoint/(this.board.stretchX);
    
    var checkPoint = new JXG.Coords(JXG.COORDS_BY_SCREEN, [x,y], this.board);
    var r = this.Radius();
    
    var dist = Math.sqrt(Math.pow(this.midpoint.coords.usrCoords[1]-checkPoint.usrCoords[1],2) + 
                         Math.pow(this.midpoint.coords.usrCoords[2]-checkPoint.usrCoords[2],2));
   
    var has = (dist < r);
    if(has) {
        var p = {};
        p.coords = new JXG.Coords(JXG.COORDS_BY_USER, 
                              [this.midpoint.coords.usrCoords[1], 
                               this.board.origin.usrCoords[2]/(this.board.stretchY)],
                              this.board);
        var angle1 = this.board.algebra.trueAngle(this.point2, this.midpoint, p);
        var angle2 = this.board.algebra.trueAngle(this.point3, this.midpoint, p);

        var xy = {};
        xy.coords = checkPoint;
        var angle3 = this.board.algebra.trueAngle(xy, this.midpoint, p); 
        if(angle1 >= angle2) {
            if(angle1 < angle3 || angle3 < angle2) {
                has = false;
            }
        }
        else {
            if(angle3 > angle1) {
                if(angle3 < angle2) {
                    has = false;
                }
            }
        }
    }
    return has;    
};

/**
 * Calculates the arcs radius.
 * @type float
 * @return The arcs radius
 */
JXG.Arc.prototype.Radius = function() {
    return(Math.sqrt(Math.pow(this.midpoint.coords.usrCoords[1]-this.point2.coords.usrCoords[1],2) + Math.pow(this.midpoint.coords.usrCoords[2]-this.point2.coords.usrCoords[2],2)));
};

/**
  * @deprecated
  */
JXG.Arc.prototype.getRadius = function() {
    this.Radius();
};

/**
 * return TextAnchor
 */
JXG.Arc.prototype.getTextAnchor = function() {
    return this.midpoint.coords;
};

/**
 * return LabelAnchor
 */
JXG.Arc.prototype.getLabelAnchor = function() {
    var angle = this.board.algebra.trueAngle(this.point2, this.midpoint, this.point3);
    var dx = 10/(this.board.stretchX);
    var dy = 10/(this.board.stretchY);
    
    var bxminusax = this.point2.coords.usrCoords[1] - this.midpoint.coords.usrCoords[1];
    var byminusay = this.point2.coords.usrCoords[2] - this.midpoint.coords.usrCoords[2];

    if(this.label.content != null) {                          
        this.label.content.relativeCoords = new JXG.Coords(JXG.COORDS_BY_USER, [0/(this.board.stretchX),0/(this.board.stretchY)],this.board);                      
    }  

    var coords = new JXG.Coords(JXG.COORDS_BY_USER, 
                          [this.midpoint.coords.usrCoords[1]+ Math.cos(angle*Math.PI/(2*180))*bxminusax - Math.sin(angle*Math.PI/(2*180))*byminusay, 
                           this.midpoint.coords.usrCoords[2]+ Math.sin(angle*Math.PI/(2*180))*bxminusax + Math.cos(angle*Math.PI/(2*180))*byminusay], 
                          this.board);

    var vecx = coords.usrCoords[1] - this.midpoint.coords.usrCoords[1];
    var vecy = coords.usrCoords[2] - this.midpoint.coords.usrCoords[2];
    
    var length = Math.sqrt(vecx*vecx+vecy*vecy);
    vecx = vecx*(length+dx)/length;
    vecy = vecy*(length+dy)/length;

    var coords2 = new JXG.Coords(JXG.COORDS_BY_USER, [this.midpoint.coords.usrCoords[1]+vecx,this.midpoint.coords.usrCoords[2]+vecy],this.board);
    
    return coords2;
};

/**
 * Uses the boards renderer to update the arc.
 * update() is not needed for arc.
 */
JXG.Arc.prototype.updateRenderer = function () {
    if (this.needsUpdate) { 
        this.board.renderer.updateArc(this);
        this.needsUpdate = false;
    }
    
    /* Update the label if visible. */
    if(this.hasLabel && this.label.content.visProp['visible'] && this.isReal) {
        //this.label.setCoordinates(this.coords);
        this.label.content.update();
        //this.board.renderer.updateLabel(this.label);
        this.board.renderer.updateText(this.label.content);
    }      
};

/**
 * Determines whether the arc has arrows at start or end of the arc.
 * @param {bool} firstArrow True if there is an arrow at the start of the arc, false otherwise.
 * @param {bool} lastArrow True if there is an arrow at the end of the arc, false otherwise.
 * Is stored at visProp['firstArrow'] and visProp['lastArrow']
 */
JXG.Arc.prototype.setArrow = function (firstArrow, lastArrow) {
    this.visProp['firstArrow'] = firstArrow;
    this.visProp['lastArrow'] = lastArrow;
     
    this.board.renderer.updateArc(this);
    
    if(this.hasLabel && this.label.content.visProp['visible']) {
        //this.label.setCoordinates(this.coords);
        this.label.content.update();
        //this.board.renderer.updateLabel(this.label);
        this.board.renderer.updateText(this.label.content);
    }     
};

/**
 * Creates a new arc.
 * @param {JXG.Board} board The board the arc is put on.
 * @param {Array} parents Array of three points defining the arc.
 * @param {Object} attributs Object containing properties for the element such as stroke-color and visibility. See @see JXG.GeometryElement#setProperty
 * @type JXG.Arc
 * @return Reference to the created arc object.
 */
JXG.createArc = function(board, parents, attributes) {
    var el;
    
    attributes = JXG.checkAttributes(attributes,{withLabel:JXG.readOption(board.options,'arc','withLabel'), layer:null});
    // Alles 3 Punkte?
    if ( (JXG.isPoint(parents[0])) && (JXG.isPoint(parents[1])) && (JXG.isPoint(parents[2]))) {
        el = new JXG.Arc(board, parents[0], parents[1], parents[2], attributes['id'], attributes['name'],attributes['withLabel'],attributes['layer']);
    } // Ansonsten eine fette Exception um die Ohren hauen
    else
        throw new Error("JSXGraph: Can't create Arc with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "' and '" + (typeof parents[2]) + "'.");

    return el;
};

JXG.JSXGraph.registerElement('arc', JXG.createArc);

/**
 * Creates a new semicircle. The semicircle is drawn clock-wise between the first and the second defining point.
 * @param {JXG.Board} board The board the semicircle is put on.
 * @param {Array} parents Array of two opposite points defining the semicircle.
 * @param {Object} attributs Object containing properties for the element such as stroke-color and visibility. See @see JXG.GeometryElement#setProperty
 * @type JXG.Arc
 * @return Reference to the created arc object.
 */
JXG.createSemicircle = function(board, parents, attributes) {
    var el, mp, idmp;
    
    attributes = JXG.checkAttributes(attributes,{withLabel:JXG.readOption(board.options,'arc','withLabel'), layer:null});
    if(attributes['id'] != null) {
        idmp = attributes['id']+'_mp';
    }
    // Alles 2 Punkte?
    if ( (JXG.isPoint(parents[0])) && (JXG.isPoint(parents[1])) ) {
        mp = board.create('midpoint', [parents[0], parents[1]], {id:idmp, withLabel:false, visible:false});
        el = new JXG.Arc(board, mp, parents[1], parents[0], attributes['id'], attributes['name'],attributes['withLabel'],attributes['layer']);
    } // Ansonsten eine fette Exception um die Ohren hauen
    else
        throw new Error("JSXGraph: Can't create Semicircle with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");

    return el;
};

JXG.JSXGraph.registerElement('semicircle', JXG.createSemicircle);

/**
 * Creates a new circumcircle arc through three defining points.
 * @param {JXG.Board} board The board the arc is put on.
 * @param {Array} parents Array of three points defining the circumcircle arc.
 * @param {Object} attributs Object containing properties for the element such as stroke-color and visibility. See @see JXG.GeometryElement#setProperty
 * @type JXG.Arc
 * @return Reference to the created arc object.
 */
JXG.createCircumcircleArc = function(board, parents, attributes) {
    var el, mp, idmp, det;
    
    attributes = JXG.checkAttributes(attributes,{withLabel:JXG.readOption(board.options,'arc','withLabel'), layer:null});
    if(attributes['id'] != null) {
        idmp = attributes['id']+'_mp';
    }
    
    // Alles 3 Punkte?
    if ( (JXG.isPoint(parents[0])) && (JXG.isPoint(parents[1])) && (JXG.isPoint(parents[2]))) {
        mp = board.create('circumcirclemidpoint',[parents[0], parents[1], parents[2]], {id:idmp, withLabel:false, visible:false});
        det = (parents[0].coords.usrCoords[1]-parents[2].coords.usrCoords[1])*(parents[0].coords.usrCoords[2]-parents[1].coords.usrCoords[2]) -
              (parents[0].coords.usrCoords[2]-parents[2].coords.usrCoords[2])*(parents[0].coords.usrCoords[1]-parents[1].coords.usrCoords[1]);
        if(det < 0) {
            el = new JXG.Arc(board, mp, parents[0], parents[2], attributes['id'], attributes['name'],attributes['withLabel'],attributes['layer']);
        }
        else {
            el = new JXG.Arc(board, mp, parents[2], parents[0], attributes['id'], attributes['name'],attributes['withLabel'],attributes['layer']);         
        }
        
        el.update = function() {
            var determinante;
            if(this.traced) {
                this.cloneToBackground(true);
            }
            determinante = (parents[0].coords.usrCoords[1]-parents[2].coords.usrCoords[1])*(parents[0].coords.usrCoords[2]-parents[1].coords.usrCoords[2]) -
                           (parents[0].coords.usrCoords[2]-parents[2].coords.usrCoords[2])*(parents[0].coords.usrCoords[1]-parents[1].coords.usrCoords[1]);
            if(determinante < 0) {
                this.point2 = parents[0];
                this.point3 = parents[2];
            }
            else {
                this.point2 = parents[2];
                this.point3 = parents[0];
            }
        };
    } // Ansonsten eine fette Exception um die Ohren hauen
    else
        throw new Error("JSXGraph: create Circumcircle Arc with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "' and '" + (typeof parents[2]) + "'.");


    return el;
};

JXG.JSXGraph.registerElement('circumcirclearc', JXG.createCircumcircleArc);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.

*/
/**
 * Creates a new instance of Sector.
 * @class Sector stores all style and functional properties that are required
 * to draw a sector on a board.
 * @param {JXG.Board} board Reference to the board the sector is drawn on.
 * @param {JXG.Point} p1 Midpoint of the sector.
 * @param {JXG.Point} p2 Point defining the sectors radius
 * @param {JXG.Point} p3 This point defines the angle of the sectors section.
 * @param {Array} ids Unique identifiers for the derived objects (arc, endpoint of the arc, line from p1 to p2, line from p1 to the endpoint of the arc) . 
 * Must be Strings. If null or an empty string is given to any of these, an unique id will be generated.
 * @param {Array} names Names for the derived objects (arc, endpoint of the arc, line from p1 to p2, line from p1 to the endpoint of the arc) . 
 * Must be Strings. If null or an empty string is given to any of these, an unique id will be generated.
 * @param {String} id Unique identifier for this object.  If null or an empty string is given,
 * an unique id will be generated by Board
 * @param {String} name Not necessarily unique name, displayed on the board.  If null or an
 * empty string is given, an unique name will be generated.
 * @see JXG.Board#addSector
 * @constructor
 * @extends JXG.GeometryElement
 */

 /* Sector legt nur die benoetigten Unterelemente an und verwaltet diese als Kinder, wird nicht mehr direkt gezeichnet */
JXG.Sector = function (board, p1, p2, p3, ids, names, id, layer) {
    /* Call the constructor of GeometryElement */
    this.constructor();
    /**
     * Sets type of GeometryElement, value is OBJECT_TYPE_SECTOR.
     * @final
     * @type int
     */    
    this.type = JXG.OBJECT_TYPE_SECTOR;
    this.elementClass = JXG.OBJECT_CLASS_AREA;                

    this.init(board, id, '');
    /**
     * Set the display layer.
     */
    if (layer == null) layer = board.options.layer['sector'];
    this.layer = layer;
    
    if(!JXG.isArray(ids)) {
        ids = [null, null, null, null];
    }
    
    if(!JXG.isArray(names)) {
     //   names = [null, null, null, null];
    }

    /**
     * Midpoint of the sector.
     * @type JXG.Point
     */
    this.point1 = JXG.getReference(this.board, p1);

    /**
     * Point defining the sectors circle.
     * @type JXG.Point
     */    
    this.point2 = JXG.getReference(this.board, p2);
    
    /**
     * The point defining the angle of the sector.
     * @type JXG.Point
     */
    this.point3 = JXG.getReference(this.board, p3);

    /**
     * This is just for the hasPoint() method. Precision for highlighting.
     * @type int
     */    
    //this.r = this.board.options.precision.hasPoint;
  
    this.visProp['visible'] = true;

    var circle = {}; // um projectToCircle benutzen zu koennen
    circle.midpoint = this.point1;
    var radius = this.Radius();
    circle.Radius = function() {
        return radius;
    };
    //-----------------
    // deprecated:
    circle.getRadius = function() {
        return radius;
    };
    //-----------------
    var p4coords = this.board.algebra.projectPointToCircle(this.point3,circle);
    
    var p = new JXG.Point(board, [p4coords.usrCoords[1], p4coords.usrCoords[2]], ids[1], names[1], true);
    p.fixed = true;
    this.addChild(p);
    p.update = function() {
        var circle = {}; // um projectToCircle benutzen zu koennen
        circle.midpoint = JXG.getReference(this.board, p1);
        var radius = (Math.sqrt(Math.pow(JXG.getReference(this.board, p1).coords.usrCoords[1]-JXG.getReference(this.board, p2).coords.usrCoords[1],2) + Math.pow(JXG.getReference(this.board, p1).coords.usrCoords[2]-JXG.getReference(this.board, p2).coords.usrCoords[2],2)));

        circle.Radius = function() {
            return radius;
        };
        //-------------------
        // deprecated
        circle.getRadius = function() {
            return radius;
        };
        //-------------------
        p4coords = this.board.algebra.projectPointToCircle(JXG.getReference(this.board, p3),circle);
        this.coords = p4coords;
        this.board.renderer.updatePoint(this);
        
        // Label mitschieben
        if(this.label.content.visProp['visible']) {
            //this.label.setCoordinates(this.coords);
            //this.board.renderer.updateLabel(this.label);
            this.label.content.update();
        }
        
        //for(var Element in this.childElements) {
        //    this.childElements[Element].update();
        //}     
    };
   
    var l1 = new JXG.Line(board, p1, p2, ids[2], names[2]);
    var l2 = new JXG.Line(board, p1, p.id, ids[3], names[3]);
    l1.setStraight(false,false);
    l2.setStraight(false,false);
   
    var a = new JXG.Arc(board, p1, p2, p3, ids[0], names[0]);
    a.visProp['fillColor'] = this.board.options.sector.fillColor;
    a.visProp['highlightFillColor'] = this.board.options.sector.highlightFillColor;
    a.visProp['fillOpacity'] = this.board.options.sector.fillOpacity;
    a.visProp['highlightFillOpacity'] = this.board.options.sector.highlightFillOpacity;
    
    /**
     * Endpoint of the derived arc.
     * @type JXG.Point
     */
    this.point4 = p;
    
    /**
     * The derived lines. 
    * this.lines[0] is the line between the midpoint of the sector and the startpoint of the arc
    * this.lines[1] is the line between the midpoint of the sector and the endpoint of the arc
     * @type Array
     */    
    this.lines = [l1, l2];

    /**
     * The derived arc.
     * @type JXG.Arc
     */    
    this.arc = a;
    
    /* Register sector at board */
    this.id = this.board.addSector(this);
    
    /* Add sector as child to defining points */
    this.point1.addChild(this);
    this.point2.addChild(this);
    this.point3.addChild(this);
    
    return this;
};   
JXG.Sector.prototype = new JXG.GeometryElement;

/**
 * Checks whether (x,y) is near the sector.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @return {bool} Always false, because the sectors interior shall not be highlighted
 */
JXG.Sector.prototype.hasPoint = function (x, y) { 
    return false; 
};

/**
 * Calculates the sectors radius.
 * @type float
 * @return The sectors radius
 */
JXG.Sector.prototype.Radius = function() {
    return(Math.sqrt(Math.pow(this.point1.coords.usrCoords[1]-this.point2.coords.usrCoords[1],2) + Math.pow(this.point1.coords.usrCoords[2]-this.point2.coords.usrCoords[2],2)));
};

/**
 *@deprecated
 */
JXG.Sector.prototype.getRadius = function() {
    return this.Radius();
};

/**
 * Uses the boards renderer to update the sector and all of its children.
 */
 JXG.Sector.prototype.updateRenderer = function () {
   /* nichts zu tun */
};

JXG.createSector = function(board, parentArr, atts) {
    var el;
    atts = JXG.checkAttributes(atts,{withLabel:JXG.readOption(board.options,'sector','withLabel'), layer:null});
    // Alles 3 Punkte?
    if ( (JXG.isPoint(parentArr[0])) && (JXG.isPoint(parentArr[1])) && (JXG.isPoint(parentArr[2]))) {
        el = new JXG.Sector(board, parentArr[0], parentArr[1], parentArr[2], atts["ids"], atts["names"], atts['id'], atts['layer']);
    } // Ansonsten eine fette Exception um die Ohren hauen
    else
        throw new Error("JSXGraph: Can't create sector with parent types '" + (typeof parentArr[0]) + "' and '" + (typeof parentArr[1]) + "' and '" + (typeof parentArr[2]) + "'.");

    return el;
};

JXG.JSXGraph.registerElement('sector', JXG.createSector);

/**
 * Creates a new circumcircle sector through three defining points.
 * @param {JXG.Board} board The board the sector is put on.
 * @param {Array} parents Array of three points defining the circumcircle sector.
 * @param {Object} attributs Object containing properties for the element such as stroke-color and visibility. See @see JXG.GeometryElement#setProperty
 * @type JXG.Sector
 * @return Reference to the created arc object.
 */
 JXG.createCircumcircleSector = function(board, parents, attributes) {
    var el, mp, idmp, det;
    
    attributes = JXG.checkAttributes(attributes,{withLabel:JXG.readOption(board.options,'sector','withLabel'), layer:null});
    if(attributes['id'] != null) {
        idmp = attributes['id']+'_mp';
    }
    
    // Alles 3 Punkte?
    if ( (JXG.isPoint(parents[0])) && (JXG.isPoint(parents[1])) && (JXG.isPoint(parents[2]))) {
        mp = board.create('circumcirclemidpoint',[parents[0], parents[1], parents[2]], {id:idmp, withLabel:false, visible:false});
        det = (parents[0].coords.usrCoords[1]-parents[2].coords.usrCoords[1])*(parents[0].coords.usrCoords[2]-parents[1].coords.usrCoords[2]) -
              (parents[0].coords.usrCoords[2]-parents[2].coords.usrCoords[2])*(parents[0].coords.usrCoords[1]-parents[1].coords.usrCoords[1]);
        if(det < 0) {
            el = new JXG.Sector(board, mp, parents[0], parents[2], attributes['id'], [attributes['name'],'','',''],attributes['withLabel'],attributes['layer']);
        }
        else {
            el = new JXG.Sector(board, mp, parents[2], parents[0], attributes['id'], [attributes['name'],'','',''],attributes['withLabel'],attributes['layer']);         
        }
        
        el.arc.update = function() {
            var determinante;
            if(this.traced) {
                this.cloneToBackground(true);
            }
            determinante = (parents[0].coords.usrCoords[1]-parents[2].coords.usrCoords[1])*(parents[0].coords.usrCoords[2]-parents[1].coords.usrCoords[2]) -
                           (parents[0].coords.usrCoords[2]-parents[2].coords.usrCoords[2])*(parents[0].coords.usrCoords[1]-parents[1].coords.usrCoords[1]);
            if(determinante < 0) {
                this.point2 = parents[0];
                this.point3 = parents[2];
            }
            else {
                this.point2 = parents[2];
                this.point3 = parents[0];
            }    
        }
        el.point4.setProperty({visible:false});
    } // Ansonsten eine fette Exception um die Ohren hauen
    else
        throw new Error("JSXGraph: Can't create circumcircle sector with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "' and '" + (typeof parents[2]) + "'.");

    return el;
};

JXG.JSXGraph.registerElement('circumcirclesector', JXG.createCircumcircleSector);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/** 
 * @fileoverview This file contains the class definition of JXG.Angle and the JXG.createAngle
 * element wrapper.
 * @author graphjs
 */
 
/**
 * Creates a new instance of Angle.
 * @class JXG.Angle holds properties and methods for creating and modifying the geometry
 * element angle which is visualized with a sector.
 * @param {JXG.Board} board The board the angle is associated with.
 * @param {JXG.Point} p1 First point A defining the angle ABC.
 * @param {JXG.Point} p2 Second point B defining the angle ABC.
 * @param {JXG.Point} p3 Third point C defining the angle ABC.
 * @param {number} radius Radius of the angle.
 * @param {String} text A text drawn beside the angle.
 * @param {String} id Unique identifier for this object. If null or an empty string is given,
 * an unique id will be generated by Board
 * @param {String} name Not necessarily unique name, displayed on the board. If null or an
 * empty string is given, an unique name will be generated.
 * @constructor
 * @extends JXG.GeometryElement
 */
JXG.Angle = function (board, p1, p2, p3, radius, text, id, name, withLabel, layer) {
    /* Call the constructor of GeometryElement */
    this.constructor();
    /**
     * Type of GeometryElement, value is OBJECT_TYPE_ANGLE.
     * @constant
     * @type number
     */    
    this.type = JXG.OBJECT_TYPE_ANGLE;
    
    /**
     * Class of the element, value is OBJECT_CLASS_AREA.
     * @constant
     * @type number
     */
    this.elementClass = JXG.OBJECT_CLASS_AREA;
    
    this.init(board, id, name);

    /**
     * Set the display layer.
     */
    if (layer == null) layer = board.options.layer['angle'];
    this.layer = layer;

    /**
     * First point A defining the angle ABC. Do no set this property directly as it
     * will break JSXGraph's dependency tree.
     * @type JXG.Point
     * @private
     */
    this.point1 = JXG.getReference(this.board, p1);

    /**
     * Second point B defining the angle ABC. Do no set this property directly as it
     * will break JSXGraph's dependency tree.
     * @type JXG.Point
     * @private
     */    
    this.point2 = JXG.getReference(this.board, p2);
    
    /**
     * Third point C defining the angle ABC. Do no set this property directly as it
     * will break JSXGraph's dependency tree.
     * @type JXG.Point
     * @private
     */
    this.point3 = JXG.getReference(this.board, p3);    

    /**
    * Determines the radius of the sector that visualizes the angle.
    * @type number
    */
    this.radius = this.board.options.angle.radius;
    if(radius != undefined && radius != null) {
        this.radius = radius;
    }

    this.visProp['fillColor'] = this.board.options.angle.fillColor;
    this.visProp['highlightFillColor'] = this.board.options.angle.highlightFillColor;
    this.visProp['fillOpacity'] = this.board.options.angle.fillOpacity;
    this.visProp['highlightFillOpacity'] = this.board.options.angle.highlightFillOpacity;
    this.visProp['strokeColor'] = this.board.options.angle.strokeColor;    


    /* TODO Why is this here and not in board.generateName? --michael */
    if(text == '') {
        var possibleNames = ['&alpha;', '&beta;', '&gamma;', '&delta;', '&epsilon;', '&zeta;', '&eta', '&theta;',
                                '&iota;', '&kappa;', '&lambda;', '&mu;', '&nu;', '&xi;', '&omicron;', '&pi;', '&rho;', 
                                '&sigmaf;', '&sigma;', '&tau;', '&upsilon;', '&phi;', '&chi;', '&psi;', '&omega;'],
            i = 0,
            j, x, el, pre, post, found;
        while(i < possibleNames.length) {
            j=i;
            x = possibleNames[i];
            for(el in board.objects) {
                if(board.objects[el].type == JXG.OBJECT_TYPE_ANGLE) {
                    if(board.objects[el].text == x) {
                        i++;
                        break;
                    }
                }
            }
            if(i==j) {
                text = x;
                i = possibleNames.length+1;
            }
        }
        if(i == possibleNames.length) {
            pre = '&alpha;_{';
            post = '}';
            found = false;
            j=0;
            while(!found) {
                for(el in board.objects) {
                    if(board.objects[el].type == JXG.OBJECT_TYPE_ANGLE) {
                        if(board.objects[el].text == (pre+j+post)) {
                            found = true;
                            break;
                        }
                    }
                }
                if(found) {
                    found= false;
                }
                else {
                    found = true;
                    text = (pre+j+post);
                }
            }
        }
    }

    /** 
    * Text (i.e. name) of the Angle.
    * @type String
    * @private
    */    
    this.text = text;

    // create Label
    var tmp = this.name;
    this.name = this.text;
    this.createLabel(withLabel);    
    this.name = tmp;
    
    this.id = this.board.addAngle(this);
    
    /* Add sector as child to defining points */
    this.point1.addChild(this);
    this.point2.addChild(this);
    this.point3.addChild(this);    
};

JXG.Angle.prototype = new JXG.GeometryElement;

/**
 * Checks whether (x,y) is near the angle. This method is a stub always returning
 * false.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @return {bool} Always false, because the angles interior shall not be highlighted
 */
JXG.Angle.prototype.hasPoint = function (x, y) { 
    return false; 
};

/**
 * Uses the boards renderer to update the angle and all of its children.
 */
 JXG.Angle.prototype.updateRenderer = function () {
    if (this.needsUpdate) {
        this.board.renderer.updateAngle(this);
        this.needsUpdate = false;
    }
    
    /* Update the label if visible. */
    if(this.hasLabel && this.label.content.visProp['visible'] && this.isReal) {
        //this.label.setCoordinates(this.coords);
        this.label.content.update();
        //this.board.renderer.updateLabel(this.label);
        this.board.renderer.updateText(this.label.content);
    }      
};

/**
 * return LabelAnchor
 */
JXG.Angle.prototype.getLabelAnchor = function() {
    var angle = this.board.algebra.trueAngle(this.point1, this.point2, this.point3);
    var dist = this.point1.coords.distance(JXG.COORDS_BY_USER,this.point2.coords);
    var bxminusax = (this.point1.coords.usrCoords[1] - this.point2.coords.usrCoords[1])*(this.radius/2)/dist;
    var byminusay = (this.point1.coords.usrCoords[2] - this.point2.coords.usrCoords[2])*(this.radius/2)/dist;
    var c = new JXG.Coords(JXG.COORDS_BY_USER, 
                          [this.point2.coords.usrCoords[1]+ Math.cos(angle*Math.PI/(2*160))*bxminusax - Math.sin(angle*Math.PI/(2*160))*byminusay, 
                           this.point2.coords.usrCoords[2]+ Math.sin(angle*Math.PI/(2*160))*bxminusax + Math.cos(angle*Math.PI/(2*160))*byminusay], 
                          this.board);
    if(this.label.content != null) {                          
        this.label.content.relativeCoords = new JXG.Coords(JXG.COORDS_BY_USER, [0/(this.board.stretchX),0/(this.board.stretchY)],this.board);                      
    }
    return c;
};

/**
 * Creates a new angle.
 * @param {JXG.Board} board The board the angle is put on.
 * @param {Array} parents Array of three points defining the angle.
 * @param {Object} attributs Object containing properties for the element such as stroke-color and visibility. @see JXG.GeometryElement#setProperty
 * @type JXG.Angle
 * @return Reference to the created angle object.
 */
JXG.createAngle = function(board, parents, attributes) {
    var el;
    
    attributes = JXG.checkAttributes(attributes,{withLabel:JXG.readOption(board.options,'angle','withLabel'), text:'', layer:null});
    
    // Alles 3 Punkte?
    if ( (JXG.isPoint(parents[0])) && (JXG.isPoint(parents[1])) && (JXG.isPoint(parents[2]))) {
        el = new JXG.Angle(board, parents[0], parents[1], parents[2], attributes['radius'], attributes['text'], attributes['id'], attributes['name'],attributes['withLabel'],attributes['layer']);
    } // Ansonsten eine fette Exception um die Ohren hauen
    else
        throw new Error("JSXGraph: Can't create angle with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "' and '" + (typeof parents[2]) + "'.");

    return el;
};

JXG.JSXGraph.registerElement('angle', JXG.createAngle);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph. If not, see <http://www.gnu.org/licenses/>.
*/

/** 
 * @fileoverview This file contains the class Algebra, a class for calculating algebraic/geometric
 * stuff like intersection points, angles, midpoint, and so on.
 */
 
/**
 * Creates a new instance of Algebra.
 * @class A class for algebraic computations like determining intersection points, angles, midpoints, ...
 * @param board The board the algebra object is associated with.
 * @constructor
 */
JXG.Algebra = function (/** JXG.Board */ board) {
    /**
     * Reference to board.
     * @type JXG.Board
     */
    this.board = board;
    
    /**
     * Defines float precision. Every number <tt>f</tt> with
     * Math.abs(f) < eps is assumed to be zero.
     * @default {@link JXG.Math#eps}
     * @see JXG.Math#eps
     */
    this.eps = JXG.Math.eps;
};

/**
 * Calculates the angle defined by the points A, B, C.
 * @param {JXG.Point,array} A A point  or [x,y] array.
 * @param {JXG.Point,array} B Another point or [x,y] array.
 * @param {JXG.Point,array} C A circle - no, of course the third point or [x,y] array.
 * @type number
 * @return The angle in radian measure.
 * @deprecated Use {@link JXG.Algebra#rad} instead.
 * @see #rad
 * @see #trueAngle
 */
JXG.Algebra.prototype.angle = function(A, B, C) {   
    var a = [],
        b = [],
        c = [],
        u, v, s, t;
        
    if (A.coords == null) {
        a[0] = A[0];
        a[1] = A[1];
    } else {
        a[0] = A.coords.usrCoords[1];
        a[1] = A.coords.usrCoords[2];
    }
    if (B.coords == null) {
        b[0] = B[0];
        b[1] = B[1];
    } else {
        b[0] = B.coords.usrCoords[1];
        b[1] = B.coords.usrCoords[2];
    }
    if (C.coords == null) {
        c[0] = C[0];
        c[1] = C[1];
    } else {
        c[0] = C.coords.usrCoords[1];
        c[1] = C.coords.usrCoords[2];
    }
    u = a[0] - b[0];
    v = a[1] - b[1];
    s = c[0] - b[0];
    t = c[1] - b[1];
    return Math.atan2(u*t-v*s,u*s+v*t);    
};

/**
 * Calculates the angle defined by the three points A, B, C if you're going from A to C around B counterclockwise.
 * @param A Point or [x,y] array
 * @param B Point or [x,y] array
 * @param C Point or [x,y] array
 * @return The angle in degrees.
 * @see #rad
 */
JXG.Algebra.prototype.trueAngle = function(/** JXG.Point */ A, /** JXG.Point */ B, /** JXG.Point */ C) /** number */ {
    return this.rad(A,B,C)*57.295779513082323; // *180.0/Math.PI;
};

/**
 * Calculates the internal angle defined by the three points A, B, C if you're going from A to C around B counterclockwise.
 * @param {JXG.Point} A Point or [x,y] array
 * @param {JXG.Point} B Point or [x,y] array
 * @param {JXG.Point} C Point or [x,y] array
 * @type number
 * @see #trueAngle
 * @return Angle in radians.
 */
JXG.Algebra.prototype.rad = function(A,B,C) {
    var ax, ay, bx, by, cx, cy,
        abx, aby, cbx, cby,
        cp, l1, l2, phiacos, phicos, sp, 
        phi = 0;
        
    if (A.coords == null) {
        ax = A[0];
        ay = A[1];
    } else {
        ax = A.coords.usrCoords[1];
        ay = A.coords.usrCoords[2];
    }
    if (B.coords == null) {
        bx = B[0];
        by = B[1];
    } else {
        bx = B.coords.usrCoords[1];
        by = B.coords.usrCoords[2];
    }
    if (C.coords == null) {
        cx = C[0];
        cy = C[1];
    } else {
        cx = C.coords.usrCoords[1];
        cy = C.coords.usrCoords[2];
    }
    cbx = cx - bx;
    cby = cy - by;
    abx = ax - bx;
    aby = ay - by;
    
    sp = cbx*abx + cby*aby;               // scalar product of c-b and a-b
    cp = abx*cby - aby*cbx;               // cross product of a-b c-b
    l1 = Math.sqrt(abx*abx + aby*aby);    // length of a-b
    l2 = Math.sqrt(cbx*cbx + cby*cby);    // length of c-b
    phiacos = sp / (l1 * l2);             // calculate the angle as cosine from scalar product
    if (phiacos > 1) { // these things should not happen, but can happen because of numerical inaccurracy
        phiacos = 1;
    } else if (phiacos < -1) {
        phiacos = -1;
    }
    phicos = Math.acos(phiacos); // calculate the angle
        /*
         * The calculated angle may not be the right angle because of the calculation of acos 
        real     | quadrant  | quadrant | algebraic sign 
        quadrant | by cosine | by sine  | of sine 
           1.    |   1.      |   1.     |   +
           2.    |   2.      |   1.     |   +
           3.    |   2.      |   3.     |   -
           4.    |   1.      |   3.     |   - 
         So only for the first quadrant the calculated angle is ok. 
         But we can use the sine, which is connected with the cross product to select the right angle. 
         Calculate the sine of the calculated angle and multiply it with the cross product's value.
         
        real     | quadrant  | algebraic sign | algebraic sign of 
        quadrant | by cosine | of sin(phicos) | cross product
           1.    |   1.      |   +            |   +
           2.    |   2.      |   +            |   +
           3.    |   2.      |   +            |   -
           4.    |   1.      |   +            |   - 
         So always the negative angle of phicos has to be taken if the product is negative.
         */
    if ((Math.sin(phicos) * cp) < 0) {
        phi = 6.2831853071795862 - phicos; // 2 * Math.PI - phicos;
    } else {
        phi = phicos;
    }
    return phi;
};

/**
 * Calculates the bisection between the three points A, B, C. The bisection is defined by two points:
 * Parameter B and a point with the coordinates calculated in this function.
 * @param A Point
 * @param B Point
 * @param C Point
 * @return Coordinates of the second point defining the bisection.
 */
JXG.Algebra.prototype.angleBisector = function(/** JXG.Point */ A, /** JXG.Point */ B, /** JXG.Point */ C) /** JXG.Coords */ {
    /* First point */
    var Ac = A.coords.usrCoords,
        Bc = B.coords.usrCoords, 
        Cc = C.coords.usrCoords,
        x = Ac[1]-Bc[1],
        y = Ac[2]-Bc[2],
        d = Math.sqrt(x*x+y*y),
        phiA, phiC, phi;
    x /= d;
    y /= d;
    
    phiA = Math.acos(x);
    if (y<0) { phiA *= -1; }
    if (phiA<0) { phiA += 2*Math.PI; } 
    
    /* Second point */
    x = Cc[1]-Bc[1];
    y = Cc[2]-Bc[2];
    d = Math.sqrt(x*x+y*y);
    x /= d;
    y /= d;
    
    phiC = Math.acos(x);
    if (y<0) { phiC *= -1; }
    if (phiC<0) { phiC += 2*Math.PI; } 
 
    phi=(phiA+phiC)*0.5;
    if (phiA>phiC) { 
        phi+=Math.PI;
    }

    x = Math.cos(phi)+Bc[1];
    y = Math.sin(phi)+Bc[2];
    
    return new JXG.Coords(JXG.COORDS_BY_USER, [x,y], this.board);
};
    
/**
 @private
 @deprecated
 OBSOLETE
 * Calculates the midpoint between two points A and B.
 * @param {JXG.Point} A Point
 * @param {JXG.Point} B Point
 * @type JXG.Coords
 * @return Coordinates of the point right in the middle of two given points.
 */
JXG.Algebra.prototype.midpoint = function(A, B) {   
    return new JXG.Coords(JXG.COORDS_BY_USER, 
                      [(A.coords.usrCoords[0] + B.coords.usrCoords[0])*0.5, 
                       (A.coords.usrCoords[1] + B.coords.usrCoords[1])*0.5, 
                       (A.coords.usrCoords[2] + B.coords.usrCoords[2])*0.5], 
                      this.board);
};

/**
 @private
 @deprecated
 OBSOLETE
 * Calculates the coordinates of a point on the parallel through the given point to the given line through point1 and point2.
 * @param {JXG.Point} point1 First point lying on the given line.
 * @param {JXG.Point} point2 Second point lying on the given line.
 * @param {JXG.Point} point Point through which the parallel is drawn.
 * @type JXG.Coords
 * @return Coordinates of a point defining the parallel together with the given point.
 */
JXG.Algebra.prototype.parallel = function(point1, point2, point) {
    var factor = 1,
        pc = point.coords.usrCoords,
        p1c = point1.coords.usrCoords,
        p2c = point2.coords.usrCoords,
        x = pc[1] + factor*(p2c[1] - p1c[1]),
        y = pc[2] + factor*(p2c[2] - p1c[2]);
    
    return new JXG.Coords(JXG.COORDS_BY_USER, [x,y], this.board);
};

/**
 * Reflects the point along the line.
 * @param {JXG.Line} line Axis of reflection.
 * @param {JXG.Point} point Point to reflect.
 * @type JXG.Coords
 * @return Coordinates of the reflected point.
 */  
JXG.Algebra.prototype.reflection = function(line,point) {
    /* (v,w) defines the slope of the line */    
    var pc = point.coords.usrCoords,
        p1c = line.point1.coords.usrCoords,
        p2c = line.point2.coords.usrCoords,
        x0, y0, x1, y1, v, w, mu;
        
    v = p2c[1]-p1c[1];
    w = p2c[2]-p1c[2];
    
    x0 = pc[1]-p1c[1];
    y0 = pc[2]-p1c[2];
    
    mu = (v*y0-w*x0)/(v*v+w*w);
    
    /* point + mu*(-y,x) waere Lotpunkt */
    x1 = pc[1] + 2*mu*w;
    y1 = pc[2] - 2*mu*v;
    
    return new JXG.Coords(JXG.COORDS_BY_USER, [x1,y1], this.board);
};

/**
 * Computes the new position of a point which is rotated 
 * around a second point (called rotpoint) by the angle phi.
 * @param {JXG.Point} rotpoint Center of the rotation
 * @param {JXG.Point} point point to be rotated
 * @param {number} phi rotation angle in arc length
 * @type JXG.Coords
 * @return Coordinates of the new position.
 */
JXG.Algebra.prototype.rotation = function(rotpoint, point, phi) {
    // 180 degrees:
    //var x0 = 2*rotpoint.coords.usrCoords[1]-point.coords.usrCoords[1];
    //var y0 = 2*rotpoint.coords.usrCoords[2]-point.coords.usrCoords[2];
    var pc = point.coords.usrCoords,
        rotpc = rotpoint.coords.usrCoords,
        x0, y0, c, s, x1, y1;
        
    x0 = pc[1]-rotpc[1];
    y0 = pc[2]-rotpc[2];
    
    c = Math.cos(phi);
    s = Math.sin(phi);
    
    x1 = x0*c-y0*s + rotpc[1];
    y1 = x0*s+y0*c + rotpc[2];
    
    return new JXG.Coords(JXG.COORDS_BY_USER, [x1,y1], this.board);
};

/**
 * Calculates the coordinates of a point on the perpendicular to the given line through
 * the given point.
 * @param {JXG.Line} line A line.
 * @param {JXG.Point} point Intersection point of line to perpendicular.
 * @type JXG.Coords
 * @return Coordinates of a point on the perpendicular to the given line through the given point.
 */
JXG.Algebra.prototype.perpendicular = function(line, point) {
    var A = line.point1.coords.usrCoords,
        B = line.point2.coords.usrCoords,
        C = point.coords.usrCoords,
        x, y, change,
        fmd, emc, d0, d1, den;
    
    if(point == line.point1) { // Punkt ist erster Punkt der Linie
        x = A[1] + B[2] - A[2];
        y = A[2] - B[1] + A[1];
        change = true;
    }
    else if(point == line.point2) {  // Punkt ist zweiter Punkt der Linie    
        x = B[1] + A[2] - B[2];
        y = B[2] - A[1] + B[1];
        change = false;
    }
    else if( ((Math.abs(A[1] - B[1]) > this.eps) && 
             (Math.abs(C[2] - (A[2] - B[2])*(C[1]-A[1])/(A[1] - B[1])-A[2]) < this.eps)) ||
             ((Math.abs(A[1] - B[1]) <= this.eps) && (Math.abs(A[1] - C[1]) < this.eps)) ) { // Punkt liegt auf der Linie
        x = C[1] + B[2] - C[2];
        y = C[2] - B[1] + C[1]; 
        change = true;
        if(Math.abs(x - C[1]) < this.eps && Math.abs(y - C[2]) < this.eps) {
            x = C[1] + A[2] - C[2];
            y = C[2] - A[1] + C[1];
            change = false;
        }
    }
    else { // Punkt liegt nicht auf der Linie -> als zweiter Punkt wird der Lotfusspunkt gewaehlt
        fmd = A[2] - B[2];
        emc = A[1] - B[1];
        d0 = B[1]*fmd - B[2]*emc;
        d1 = C[1]*emc + C[2]*fmd;
        den = fmd*fmd + emc*emc;
        if(Math.abs(den)<this.eps) {
            den = this.eps;
        }
        x = (d0*fmd + d1*emc) / den;
        y = (d1*fmd - d0*emc) /den;
        change = true;
    }                            
    return [new JXG.Coords(JXG.COORDS_BY_USER, [x, y], this.board),change];             
};

/**
 * Calculates the midpoint of the circumcircle of the three given points.
 * @param {JXG.Point} point1 Point
 * @param {JXG.Point} point2 Point
 * @param {JXG.Point} point3 Point
 * @type JXG.Coords
 * @return Coordinates of the midpoint of the circumcircle of the given points.
 */
JXG.Algebra.prototype.circumcenterMidpoint = function(point1, point2, point3) {
    var A = point1.coords.usrCoords,
        B = point2.coords.usrCoords,
        C = point3.coords.usrCoords,
        u, v, den, x, y;

    u = ((A[1]-B[1])*(A[1]+B[1]) + (A[2]-B[2])*(A[2]+B[2])) * 0.5;
    v = ((B[1]-C[1])*(B[1]+C[1]) + (B[2]-C[2])*(B[2]+C[2])) * 0.5;
    den = (A[1]-B[1])*(B[2]-C[2]) - (B[1]-C[1])*(A[2]-B[2]);
              
    if (Math.abs(den) < this.eps) {
        den = this.eps;
    }
    
    x = (u * (B[2]-C[2]) - v*(A[2]-B[2])) / den;
    y = (v * (A[1]-B[1]) - u*(B[1]-C[1])) / den;
    
    return new JXG.Coords(JXG.COORDS_BY_USER, [x, y], this.board);
};

/**
 * Calculates the coordinates of the intersection of the given lines.
 * @param {JXG.Line} line1 Line.
 * @param {JXG.Line} line2 Line.
 * @type JXG.Coords
 * @return Coordinates of the intersection point of the given lines.
 */
JXG.Algebra.prototype.intersectLineLine = function(line1, line2) {
    var A = line1.point1.coords.usrCoords,
        B = line1.point2.coords.usrCoords,
        C = line2.point1.coords.usrCoords,
        D = line2.point2.coords.usrCoords,
        d0, d1, den, x, y;
           
    d0 = A[1]*B[2] - A[2]*B[1];
    d1 = C[1]*D[2] - C[2]*D[1];
    den = (B[2]-A[2])*(C[1]-D[1]) - (A[1]-B[1])*(D[2]-C[2]);
                 
    if(Math.abs(den) < this.eps) {
         den = this.eps; 
    }
    x = (d0*(C[1]-D[1]) - d1*(A[1]-B[1])) / den;
    y = (d1*(B[2]-A[2]) - d0*(D[2]-C[2])) / den;

    return new JXG.Coords(JXG.COORDS_BY_USER, [x, y], this.board);
};

/**
 * Calculates the coordinates of the intersection of the given line and circle.
 * @param {JXG.Circle} circle Circle.
 * @param {JXG.Line} line Line.
 * @type array
 * @return Array of the Coordinates of the intersection points of the given circle with the given line and
 * the amount of intersection points in the first component of the array.
 */
JXG.Algebra.prototype.intersectCircleLine = function(circle, line) {
    var eA = line.point1.coords.usrCoords,
        eB = line.point2.coords.usrCoords,
        fM = circle.midpoint.coords.usrCoords,
        s, d0, d1, b, w, h, r, n1, dx, dy, firstPointX, firstPointY, l, x, y, n1s, firstPoint, secondPoint, d;

    s = line.point1.Dist(line.point2);
    if (s > 0) {
        d0 = circle.midpoint.Dist(line.point1);
        d1 = circle.midpoint.Dist(line.point2);
        b = ((d0 * d0) + (s * s) - (d1 * d1)) / (2 * s);
        w = (d0 * d0) - (b * b);
        w = (w < 0) ? 0 : w;
        h = Math.sqrt(w);
        
        r = circle.Radius();
        n1 = Math.sqrt((r * r) - h*h);
        dx = eB[1] - eA[1];
        dy = eB[2] - eA[2];
        firstPointX = fM[1] + (h / s) * dy;
        firstPointY = fM[2] - (h / s) * dx;
        d0 = (eB[1] * dy) - (eB[2] * dx);
        d1 = (firstPointX * dx) + (firstPointY * dy);
        l = (dy * dy) + (dx * dx);
        if (Math.abs(l) < this.eps) { l = this.eps; }
        x = ((d0 * dy) + (d1 * dx)) / l;
        y = ((d1 * dy) - (d0 * dx)) / l;
        n1s = n1/s;
        firstPoint =  new JXG.Coords(JXG.COORDS_BY_USER, [x + n1s * dx, y + n1s * dy], this.board);
        secondPoint = new JXG.Coords(JXG.COORDS_BY_USER, [x - n1s * dx, y - n1s * dy], this.board);
        d = circle.midpoint.coords.distance(JXG.COORDS_BY_USER, firstPoint);
      
        if ((r < (d - 1)) || isNaN(d)) {
            return [0];
        } else {
            return [2,firstPoint,secondPoint];       
        }
    }
    return [0];
};

/**
 * Calculates the coordinates of the intersection of the given circles.
 * @param {JXG.Circle} circle1 Circle.
 * @param {JXG.Circle} circle2 Circle.
 * @type array
 * @return Array of the Coordinates of the intersection points of the given circles and the
 * amount of intersection points in the first component of the array.
 */
JXG.Algebra.prototype.intersectCircleCircle = function(circle1, circle2) { 
    var intersection = {},
        r1 = circle1.Radius(),
        r2 = circle2.Radius(),
        M1 = circle1.midpoint.coords.usrCoords,
        M2 = circle2.midpoint.coords.usrCoords,
        rSum, rDiff, s, 
        dx, dy, a, h;
        
    rSum = r1 + r2;
    rDiff = Math.abs(r1 - r2);    
    // Abstand der Mittelpunkte der beiden Kreise
    s = circle1.midpoint.coords.distance(JXG.COORDS_BY_USER, circle2.midpoint.coords);
    if (s > rSum) {
        return [0]; // Kreise schneiden sich nicht, liegen nebeneinander
    } 
    else if (s < rDiff) {
        return [0]; // Kreise schneiden sich nicht, liegen ineinander
    } 
    else {
        if (s != 0) {
            intersection[0] = 1; // es gibt einen Schnitt        
            dx = M2[1] - M1[1];
            dy = M2[2] - M1[2];
            a = (s * s - r2 * r2 + r1 * r1) / (2 * s);
            h = Math.sqrt(r1 * r1 - a * a);
            intersection[1] = new JXG.Coords(JXG.COORDS_BY_USER, 
                                             [M1[1] + (a / s) * dx + (h / s) * dy, 
                                              M1[2] + (a / s) * dy - (h / s) * dx], 
                                             this.board);
            intersection[2] = new JXG.Coords(JXG.COORDS_BY_USER, 
                                             [M1[1] + (a / s) * dx - (h / s) * dy, 
                                              M1[2] + (a / s) * dy + (h / s) * dx], 
                                             this.board);    
        }
        else {
            return [0]; // vorsichtshalber... 
        }                                     
        return intersection;
    }
};

/**
 * Calculates the coordinates of the projection of a given point on a given circle. I.o.w. the
 * nearest one of the two intersection points of the line through the given point and the circles
 * midpoint.
 * @param {JXG.Point} point Point to project.
 * @param {JXG.Circle} circle Circle on that the point is projected.
 * @type JXG.Coords
 * @return The coordinates of the projection of the given point on the given circle.
 */
JXG.Algebra.prototype.projectPointToCircle = function(point,circle) {
    var dist = point.coords.distance(JXG.COORDS_BY_USER, circle.midpoint.coords),
        P = point.coords.usrCoords,
        M = circle.midpoint.coords.usrCoords,
        x, y, factor;
        
    if(Math.abs(dist) < this.eps) {
        dist = this.eps;
    }
    factor = circle.Radius() / dist;
    x = M[1] + factor*(P[1] - M[1]);
    y = M[2] + factor*(P[2] - M[2]);
    
    return new JXG.Coords(JXG.COORDS_BY_USER, [x, y], this.board);    
};

/**
 * Calculates the coordinates of the projection of a given point on a given line. I.o.w. the
 * intersection point of the given line and its perpendicular through the given point.
 * @param {JXG.Point} point Point to project.
 * @param {JXG.Line} line Line on that the point is projected.
 * @type JXG.Coords
 * @return The coordinates of the projection of the given point on the given line.
 */
JXG.Algebra.prototype.projectPointToLine = function(point, line) {
/*
    // Euclidean version
    var fmd = line.point1.coords.usrCoords[2] - line.point2.coords.usrCoords[2];
    var emc = line.point1.coords.usrCoords[1] - line.point2.coords.usrCoords[1];
    var d0 = line.point2.coords.usrCoords[1]*fmd - line.point2.coords.usrCoords[2] *emc;
    var d1 = point.coords.usrCoords[1]*emc + point.coords.usrCoords[2]*fmd;
    var den = fmd*fmd + emc*emc;
    if(Math.abs(den)<this.eps) {
        den = this.eps;
    }
    var x = (d0*fmd + d1*emc) / den;
    var y = (d1*fmd - d0*emc) /den;
    return new JXG.Coords(JXG.COORDS_BY_USER, [x,y], this.board);       
*/
    // Homogeneous version
    var v = [0,line.stdform[1],line.stdform[2]];
    v = JXG.Math.crossProduct(v,point.coords.usrCoords);
    return this.meetLineLine(v,line.stdform,0);

    //return new JXG.Coords(JXG.COORDS_BY_USER, v, this.board);       
};

/**
 * Calculates the coordinates of the projection of a given point on a given curve. 
 * Uses {@link #projectCoordsToCurve}.
 * @param {JXG.Point} point Point to project.
 * @param {JXG.Curve} graph Curve on that the point is projected.
 * @type JXG.Coords
 * @see #projectCoordsToCurve
 * @return The coordinates of the projection of the given point on the given graph.
 */
JXG.Algebra.prototype.projectPointToCurve = function(point,curve) {
    var x = point.X(),
        y = point.Y(),
        t = point.position || 0.0,
        result = this.projectCoordsToCurve(x,y,t,curve);
    point.position = result[1];      // side effect !
    return result[0];
};

/**
 * Calculates the coordinates of the projection of a coordinates pair on a given curve. In case of
 * function graphs this is the
 * intersection point of the curve and the parallel to y-axis through the given point.
 * @param {float} x coordinate to project.
 * @param {float} y coordinate to project.
 * @param {float} start value for newtons method
 * @param {JXG.Curve} graph Curve on that the point is projected.
 * @type JXG.Coords
 * @see #projectPointToCurve
 * @return Array containing the coordinates of the projection of the given point on the given graph and 
 * the position on the curve.
 */
JXG.Algebra.prototype.projectCoordsToCurve = function(x,y,t,curve) {
    var newCoords, x0, y0, x1, y1, den, i, mindist, dist, lbda,
        infty = 1000000.0;
        
    if (curve.curveType=='parameter' || curve.curveType=='polar') { 
        t = JXG.Math.Numerics.root(JXG.Math.Numerics.D(function(t){ return (x-curve.X(t))*(x-curve.X(t))+(y-curve.Y(t))*(y-curve.Y(t));}), t);
        //if (t<curve.minX()) { t = curve.minX(); }
        //if (t>curve.maxX()) { t = curve.maxX(); }
        if (t<curve.minX()) { t = curve.maxX()+t-curve.minX(); }
        if (t>curve.maxX()) { t = curve.minX()+t-curve.maxX(); }
        newCoords = new JXG.Coords(JXG.COORDS_BY_USER, [curve.X(t),curve.Y(t)], this.board);
    } else if (curve.curveType == 'plot') {
        mindist = infty;
        for (i=0;i<curve.numberPoints;i++) {
            x0 = x-curve.X(i);
            y0 = y-curve.Y(i);
            dist = Math.sqrt(x0*x0+y0*y0);
            if (dist<mindist) {
                mindist = dist;
                t = i;
            }
            if (i==curve.numberPoints-1) { continue; }

            x1 = curve.X(i+1)-curve.X(i);
            y1 = curve.Y(i+1)-curve.Y(i);
            den = x1*x1+y1*y1;
            if (den>=JXG.Math.eps) {
                lbda = (x0*x1+y0*y1)/den;
                dist = Math.sqrt( x0*x0+y0*y0 - lbda*(x0*x1+y0*y1) );
            } else {
                lbda = 0.0;
                dist = Math.sqrt(x0*x0+y0*y0);
            }
            if (lbda>=0.0 && lbda<=1.0 && dist<mindist) { 
                t = i+lbda;
                mindist = dist;
            } 
        }
        i = Math.floor(t);
        lbda = t-i;
        if (i<curve.numberPoints-1) {
            x = lbda*curve.X(i+1)+(1.0-lbda)*curve.X(i);
            y = lbda*curve.Y(i+1)+(1.0-lbda)*curve.Y(i);
        } else {
            x = curve.X(i);
            y = curve.Y(i);
        }
        newCoords = new JXG.Coords(JXG.COORDS_BY_USER, [x,y], this.board); 
    } else {             // functiongraph
        t = x;
        x = t; //curve.X(t);
        y = curve.Y(t);
        newCoords = new JXG.Coords(JXG.COORDS_BY_USER, [x,y], this.board); 
    }
    return [curve.updateTransform(newCoords),t];
};

/**
 * Calculates the coordinates of the projection of a given point on a given turtle. A turtle consists of
 * one or more curves of curveType 'plot'. Uses {@link #projectPointToCurve}.
 * @param {JXG.Point} point Point to project.
 * @param {JXG.Turtle} turtle on that the point is projected.
 * @type JXG.Coords
 * @return The coordinates of the projection of the given point on the given turtle.
 */
JXG.Algebra.prototype.projectPointToTurtle = function(point,turtle) {
    var newCoords, t, x, y, i,
        np = 0, 
        npmin = 0,
        mindist = 1000000.0, 
        dist, el, minEl, 
        len = turtle.objects.length;
    
    for(i=0;i<len;i++) {  // run through all curves of this turtle
        el = turtle.objects[i];
        if (el.type==JXG.OBJECT_TYPE_CURVE) {
            newCoords = this.projectPointToCurve(point,el);
            dist = this.distance(newCoords.usrCoords,point.coords.usrCoords);
            if (dist<mindist) {
                x = newCoords.usrCoords[1];
                y = newCoords.usrCoords[2];
                t = point.position;
                mindist = dist;
                minEl = el;
                npmin = np;
            }
            np += el.numberPoints;
        }
    }
    newCoords = new JXG.Coords(JXG.COORDS_BY_USER, [x,y], this.board);    
    point.position = t+npmin;
    return minEl.updateTransform(newCoords);
};

/**
 * Converts expression of the form <i>leftop^rightop</i> into <i>Math.pow(leftop,rightop)</i>.
 * @param {String} te Expression of the form <i>leftop^rightop</i>
 * @type String
 * @return Converted expression.
 */
JXG.Algebra.prototype.replacePow = function(te) {
    var count, pos, c,
        leftop, rightop, pre, p, left, i, right, expr;
    //te = te.replace(/\s+/g,''); // Loesche allen whitespace
                                // Achtung: koennte bei Variablennamen mit Leerzeichen
                                // zu Problemen fuehren.
    i = te.indexOf('^');
    while (i>=0) {
        left = te.slice(0,i);
        if (left.charAt(left.length-1)==')') {
            count = 1;
            pos = left.length-2;
            while (pos>=0 && count>0) {
                c = left.charAt(pos);
                if (c==')') { count++; }
                else if (c=='(') { count--; }
                pos--;
            }   
            if (count==0) {
                leftop = '';
                pre = left.substring(0,pos+1);   // finde evtl. F vor (...)^
                p = pos;
                while (p>=0 && pre.substr(p,1).match(/(\w+)/)) {
                    leftop = RegExp.$1+leftop;
                    p--;
                }
                leftop += left.substring(pos+1,left.length);
                leftop = leftop.replace(/([\(\)\+\*\%\^\-\/\]\[])/g,"\\$1");
            }
        } else {
            leftop = '\\w+';
        }
        right = te.slice(i+1);
        if (right.match(/^([\w\.]*\()/)) {
            count = 1;
            pos = RegExp.$1.length;
            while (pos<right.length && count>0) {
                c = right.charAt(pos);
                if (c==')') { count--; }
                else if (c=='(') { count++; }
                pos++;
            }
            if (count==0) {
                rightop = right.substring(0,pos);
                rightop = rightop.replace(/([\(\)\+\*\%\^\-\/\[\]])/g,"\\$1");
            }
        } else {
            rightop = '[\\w\\.]+';  // ^b 
        }
        expr = new RegExp('(' + leftop + ')\\^(' + rightop + ')');
        te = te.replace(expr,"this.board.algebra.pow($1,$2)");
        i = te.indexOf('^');
    }
    return te;
};

/**
 * Converts expression of the form <i>If(a,b,c)</i> into <i>(a)?(b):(c)/i>.
 * @param {String} te Expression of the form <i>If(a,b,c)</i>
 * @type String
 * @return Converted expression.
 */
JXG.Algebra.prototype.replaceIf = function(te) {
    var s = '',
        left, right,
        first = null,
        second = null,
        third = null,
        i, pos, count, k1, k2, c, meat;
    
    i = te.indexOf('If(');
    if (i<0) { return te; }

    te = te.replace(/""/g,'0'); // "" means not defined. Here, we replace it by 0
    while (i>=0) {
        left = te.slice(0,i);
        right = te.slice(i+3); 
        
        // Search the end of the If() command and take out the meat
        count = 1;
        pos = 0;
        k1 = -1;
        k2 = -1;
        while (pos<right.length && count>0) {
            c = right.charAt(pos);
            if (c==')') { 
                count--;
            } else if (c=='(') {
                count++;
            } else if (c==',' && count==1) {
                if (k1<0) { 
                    k1 = pos; // first komma
                } else {
                    k2 = pos; // second komma
                }
            }
            pos++;
        } 
        meat = right.slice(0,pos-1);
        right = right.slice(pos);
        
        // Test the two kommas
        if (k1<0) { return ''; } // , missing
        if (k2<0) { return ''; } // , missing
        
        first = meat.slice(0,k1);
        second = meat.slice(k1+1,k2);
        third = meat.slice(k2+1);
        first = this.replaceIf(first);    // Recurse
        second = this.replaceIf(second);  // Recurse
        third = this.replaceIf(third);    // Recurse

        s += left + '((' + first + ')?' + '('+second+'):('+third+'))';  
        te = right;
        first = null;
        second = null;
        i = te.indexOf('If(');
    }
    s += right;
    return s;
};

/**
 * Replace _{} by &lt;sub&gt;
 * @param {String} te String containing _{}.
 * @type String
 * @return Given string with _{} replaced by &lt;sub&gt;.
 */
JXG.Algebra.prototype.replaceSub = function(te) {
    if(te['indexOf']) {} else return te;

    var i = te.indexOf('_{'),
        j;
    while (i>=0) {
        te = te.substr(0,i)+te.substr(i).replace(/_\{/,'<sub>');
        j = te.substr(i).indexOf('}');
        if (j>=0) {
            te = te.substr(0,j)+te.substr(j).replace(/\}/,'</sub>');
        }
        i = te.indexOf('_{');
    }

    i = te.indexOf('_');
    while (i>=0) {
        te = te.substr(0,i)+te.substr(i).replace(/_(.?)/,'<sub>$1</sub>');
        i = te.indexOf('_');
    }
    return te;
};

/**
 * Replace ^{} by &lt;sup&gt;
 * @param {String} te String containing ^{}.
 * @type String
 * @return Given string with ^{} replaced by &lt;sup&gt;.
 */
JXG.Algebra.prototype.replaceSup = function(te) {
    if(te['indexOf']) {} else return te;

    var i = te.indexOf('^{'),
        j;
    while (i>=0) {
        te = te.substr(0,i)+te.substr(i).replace(/\^\{/,'<sup>');
        j = te.substr(i).indexOf('}');
        if (j>=0) {
            te = te.substr(0,j)+te.substr(j).replace(/\}/,'</sup>');
        }
        i = te.indexOf('^{');
    }

    i = te.indexOf('^');
    while (i>=0) {
        te = te.substr(0,i)+te.substr(i).replace(/\^(.?)/,'<sup>$1</sup>');
        i = te.indexOf('^');
    }

    return te;
};

/**
 * Replace an element's name in terms by an element's id.
 * @param term Term containing names of elements.
 * @return The same string with names replaced by ids.
 **/
JXG.Algebra.prototype.replaceNameById = function(/** string */ term) /** string */ {
    var pos = 0, end, elName, el, i,
        funcs = ['X','Y','L','V'];
    
    for (i=0;i<funcs.length;i++) {
        pos = term.indexOf(funcs[i]+'(');
        while (pos>=0) {
            if (pos>=0) {
                end = term.indexOf(')',pos+2);
                if (end>=0) {
                    elName = term.slice(pos+2,end);
                    elName = elName.replace(/\\(['"])?/g,"$1");
                    el = this.board.elementsByName[elName];
                    term = term.slice(0,pos+2) + el.id +  term.slice(end);
                }
            }
            end = term.indexOf(')',pos+2);
            pos = term.indexOf(funcs[i]+'(',end);
        }
    }

    pos = term.indexOf('Dist(');
    while (pos>=0) {
        if (pos>=0) {
            end = term.indexOf(',',pos+5);
            if (end>=0) {
                elName = term.slice(pos+5,end);
                elName = elName.replace(/\\(['"])?/g,"$1");
                el = this.board.elementsByName[elName];
                term = term.slice(0,pos+5) + el.id +  term.slice(end);
            }
        }
        end = term.indexOf(',',pos+5);
        pos = term.indexOf(',',end);
        end = term.indexOf(')',pos+1);
        if (end>=0) {
            elName = term.slice(pos+1,end);
            elName = elName.replace(/\\(['"])?/g,"$1");
            el = this.board.elementsByName[elName];
            term = term.slice(0,pos+1) + el.id +  term.slice(end);
        }
        end = term.indexOf(')',pos+1);
        pos = term.indexOf('Dist(',end);
    }

    funcs = ['Deg','Rad'];
    for (i=0;i<funcs.length;i++) {
        pos = term.indexOf(funcs[i]+'(');
        while (pos>=0) {
            if (pos>=0) {
                end = term.indexOf(',',pos+4);
                if (end>=0) {
                    elName = term.slice(pos+4,end);
                    elName = elName.replace(/\\(['"])?/g,"$1");
                    el = this.board.elementsByName[elName];
                    term = term.slice(0,pos+4) + el.id +  term.slice(end);
                }
            }
            end = term.indexOf(',',pos+4);
            pos = term.indexOf(',',end);
            end = term.indexOf(',',pos+1);
            if (end>=0) {
                elName = term.slice(pos+1,end);
                elName = elName.replace(/\\(['"])?/g,"$1");
                el = this.board.elementsByName[elName];
                term = term.slice(0,pos+1) + el.id +  term.slice(end);
            }
            end = term.indexOf(',',pos+1);
            pos = term.indexOf(',',end);
            end = term.indexOf(')',pos+1);
            if (end>=0) {
                elName = term.slice(pos+1,end);
                elName = elName.replace(/\\(['"])?/g,"$1");
                el = this.board.elementsByName[elName];
                term = term.slice(0,pos+1) + el.id +  term.slice(end);
            }
            end = term.indexOf(')',pos+1);
            pos = term.indexOf(funcs[i]+'(',end);
        }
    }
    return term;
};

/**
 * Replaces element ids in terms by element this.board.objects['id'].
 * @param term A GEONE<sub>x</sub>T function string with JSXGraph ids in it.
 * @return The input string with element ids replaced by this.board.objects["id"]. 
 **/
JXG.Algebra.prototype.replaceIdByObj = function(/** string */ term) /** string */ {
    var expr = /(X|Y|L)\(([\w_]+)\)/g;  // Suche "X(gi23)" oder "Y(gi23A)" und wandle in objects['gi23'].X() um.
    term = term.replace(expr,"this.board.objects[\"$2\"].$1()");
    
    expr = /(V)\(([\w_]+)\)/g;  // Suche "X(gi23)" oder "Y(gi23A)" und wandle in objects['gi23'].X() um.
    term = term.replace(expr,"this.board.objects[\"$2\"].Value()");

    expr = /(Dist)\(([\w_]+),([\w_]+)\)/g;  // 
    term = term.replace(expr,'this.board.objects[\"$2\"].Dist(this.board.objects[\"$3\"])');

    expr = /(Deg)\(([\w_]+),([ \w\[\w_]+),([\w_]+)\)/g;  // 
    term = term.replace(expr,'this.board.algebra.trueAngle(this.board.objects[\"$2\"],this.board.objects[\"$3\"],this.board.objects[\"$4\"])');

    expr = /Rad\(([\w_]+),([\w_]+),([\w_]+)\)/g;  // Suche Rad('gi23','gi24','gi25')
    term = term.replace(expr,'this.board.algebra.rad(this.board.objects[\"$1\"],this.board.objects[\"$2\"],this.board.objects[\"$3\"])');
    return term;
};

/**
 * Converts the given algebraic expression in GEONE<sub>x</sub>T syntax into an equivalent expression in JavaScript syntax.
 * @param {String} term Expression in GEONExT syntax
 * @type String
 * @return Given expression translated to JavaScript.
 */
JXG.Algebra.prototype.geonext2JS = function(term) {
    var expr, newterm, i,
        from = ['Abs', 'ACos', 'ASin', 'ATan','Ceil','Cos','Exp','Floor','Log','Max','Min','Random','Round','Sin','Sqrt','Tan','Trunc'], 
        to =   ['Math.abs', 'Math.acos', 'Math.asin', 'Math.atan', 'Math.ceil', 'Math.cos', 'Math.exp', 'Math.floor', 'Math.log', 'Math.max', 'Math.min', 'Math.random', 'this.board.round', 'Math.sin', 'Math.sqrt', 'Math.tan', 'Math.ceil'];
    // removed: 'Pow'  -> Math.pow
    
    //term = JXG.unescapeHTML(term);  // This replaces &gt; by >, &lt; by < and &amp; by &.ist aber zu allgemein
    term = term.replace(/&lt;/g,'<'); // Hacks, to enable not well formed XML, @see GeonextReader#replaceLessThan
    term = term.replace(/&gt;/g,'>'); 
    term = term.replace(/&amp;/g,'&'); 
    
    // Umwandeln der GEONExT-Syntax in JavaScript-Syntax
    newterm = term;
    newterm = this.replaceNameById(newterm);
    newterm = this.replaceIf(newterm);
    // Exponentiations-Problem x^y -> Math(exp(x,y).
    newterm = this.replacePow(newterm);
    newterm = this.replaceIdByObj(newterm);
    
    for (i=0; i<from.length; i++) {
        expr = new RegExp(from[i],"ig");
        newterm = newterm.replace(expr,to[i]);
    }    

    newterm = newterm.replace(/True/g,'true');
    newterm = newterm.replace(/False/g,'false');
    newterm = newterm.replace(/fasle/g,'false');

    newterm = newterm.replace(/Pi/g,'Math.PI');
    return newterm;
};

/**
 * Finds dependencies in a given term and resolves them by adding the
 * dependent object to the found objects child elements.
 * @param {JXG.GeometryElement} me Object depending on objects in given term.
 * @param {String} term String containing dependencies for the given object.
 */
JXG.Algebra.prototype.findDependencies = function(me, term) {
    var elements = this.board.elementsByName,
        el, expr, elmask;
        
    for (el in elements) {
        if (el != me.name) {
            if(elements[el].type == JXG.OBJECT_TYPE_TEXT) {
                if(!elements[el].isLabel) {
                    elmask = el.replace(/\[/g,'\\[');
                    elmask = elmask.replace(/\]/g,'\\]');
                    expr = new RegExp("\\(\(\[\\w\\[\\]'_ \]+,\)*\("+elmask+"\)\(,\[\\w\\[\\]'_ \]+\)*\\)","g");  // Searches (A), (A,B),(A,B,C)
                    if (term.search(expr)>=0) {
                        elements[el].addChild(me);
                    }
                }
            }
            else {
                elmask = el.replace(/\[/g,'\\[');
                elmask = elmask.replace(/\]/g,'\\]');
                expr = new RegExp("\\(\(\[\\w\\[\\]'_ \]+,\)*\("+elmask+"\)\(,\[\\w\\[\\]'_ \]+\)*\\)","g");  // Searches (A), (A,B),(A,B,C)
                if (term.search(expr)>=0) {
                    elements[el].addChild(me);
                }
            }
        }
    }
};

/**
 * Calculates euclidean norm for two given arrays of the same length.
 * @param {array} array1 Array of float or integer.
 * @param {array} array2 Array of float or integer.
 * @type number
 * @return Euclidean distance of the given vectors.
 */
JXG.Algebra.prototype.distance = function(array1, array2) {
    var sum = 0, 
        i, len;
        
    if(array1.length != array2.length) { return; }
    len = array1.length;
    for(i=0; i<len; i++) {
        sum += (array1[i] - array2[i])*(array1[i] - array2[i]);
    }
    return Math.sqrt(sum);
};

/**
 * Calculates euclidean distance for two given arrays of the same length.
 * If one of the arrays contains a zero in the first coordinate, and the euclidean distance
 * is different from zero it is a point at infinity and we return Infinity.
 * @param {array} array1 Array containing elements of number.
 * @param {array} array2 Array containing elements of type number.
 * @type number
 * @return Euclidean (affine) distance of the given vectors.
 */
JXG.Algebra.prototype.affineDistance = function(array1, array2) {
    var d;
    if(array1.length != array2.length) { 
        return; 
    }
    d = this.distance(array1, array2);
    if (d>this.eps && (Math.abs(array1[0])<this.eps || Math.abs(array2[0])<this.eps)) {
        return Infinity;
    } else {
        return d;
    }
};

/**
 * Compute power a^b
 * @param a Base.
 * @param b Exponent.
 * @return a to the power of b.
 */
JXG.Algebra.prototype.pow = function(/** number */ a, /** number */ b) /** number */ {
    if (a==0 || b==0) { 
        return 1;
    }
    if (Math.floor(b)==b) {// b is integer
        return Math.pow(a,b);
    } else { // b is not integer
        if (a>0) {
            return Math.exp(b*Math.log(Math.abs(a)));
        } else {
            return NaN;
        }
    }
};

/**
 * 
 * @private
 * Computes the intersection of a pair of lines, circles or both.
 * It uses the internal data array stdform of these elements.
 * @param {Array} el1 stdform of the first element (line or circle)
 * @param {Array} el2 stdform of the second element (line or circle)
 * @param {number} i Index of the intersection point that should be returned.
 * @type JXG.Coords
 * @return Coordinates of one of the possible two or more intersection points. 
 * Which point will be returned is determined by i.
 */
JXG.Algebra.prototype.meet = function(el1, el2, /** number */ i) /** JXG.Coords */ {
    var eps = this.eps; //    var eps = 0.000001;

    if (Math.abs(el1[3])<eps && Math.abs(el2[3])<eps) { // line line
        return this.meetLineLine(el1,el2,i);
    } else if (Math.abs(el1[3])>=eps && Math.abs(el2[3])<eps) { // circle line
        return this.meetLineCircle(el2,el1,i);
    } else if (Math.abs(el1[3])<eps && Math.abs(el2[3])>=eps) { // line circle
        return this.meetLineCircle(el1,el2,i);
    } else {  // circle circle
        return this.meetCircleCircle(el1,el2,i);
    }
};

/**
  * @private
  * 
  * Intersection of two lines using the stdform.
  * @param {Array} l1 stdform of the first line
  * @param {Array} l2 stdform of the second line
  * @param {number} i unused
  * @type JXG.Coords
  * @return Coordinates of the intersection point.
  */
JXG.Algebra.prototype.meetLineLine = function(l1,l2,i) {
    var s = JXG.Math.crossProduct(l1,l2);
    if (Math.abs(s[0])>this.eps) {
        s[1] /= s[0];
        s[2] /= s[0];
        s[0] = 1.0;
    }
    return new JXG.Coords(JXG.COORDS_BY_USER, s, this.board);
};

/**
  * @private
  * 
  * Intersection of line and circle using the stdform.
  * 
  * @param {Array} lin stdform of the line
  * @param {Array} circ stdform of the circle
  * @param {number} i number of the returned intersection point. 
  *   i==0: use the positive square root, 
  *   i==1: use the negative square root.
  * @type JXG.Coords
  * @return Coordinates of the intersection point
  */
 JXG.Algebra.prototype.meetLineCircle = function(lin,circ,i) {    
    var a,b,c,d,n, A,B,C, k,t;

    if (circ[4]<this.eps) { // Radius is zero, return center of circle
        return new JXG.Coords(JXG.COORDS_BY_USER, circ.slice(1,3), this.board);
    }
    c = circ[0];
    b = circ.slice(1,3);
    a = circ[3];
    d = lin[0];
    n = lin.slice(1,3);

    // Line is normalized, therefore nn==1 and we can skip some operations:
    /*
    var nn = n[0]*n[0]+n[1]*n[1];
    A = a*nn;
    B = (b[0]*n[1]-b[1]*n[0])*nn;
    C = a*d*d - (b[0]*n[0]+b[1]*n[1])*d + c*nn;
    */
    A = a;
    B = (b[0]*n[1]-b[1]*n[0]);
    C = a*d*d - (b[0]*n[0]+b[1]*n[1])*d + c;

    k = B*B-4*A*C;
    if (k>=0) {
        k = Math.sqrt(k);
        t = [(-B+k)/(2*A),(-B-k)/(2*A)];
        return ((i==0)
            ? new JXG.Coords(JXG.COORDS_BY_USER, [-t[0]*(-n[1])-d*n[0],-t[0]*n[0]-d*n[1]], this.board)
            : new JXG.Coords(JXG.COORDS_BY_USER, [-t[1]*(-n[1])-d*n[0],-t[1]*n[0]-d*n[1]], this.board)
            );
/*
            new JXG.Coords(JXG.COORDS_BY_USER, [-t[0]*(-n[1])-d*n[0]/nn,-t[0]*n[0]-d*n[1]/nn], this.board),
            new JXG.Coords(JXG.COORDS_BY_USER, [-t[1]*(-n[1])-d*n[0]/nn,-t[1]*n[0]-d*n[1]/nn], this.board)
*/
    } else {
        return new JXG.Coords(JXG.COORDS_BY_USER, [NaN,NaN], this.board);
    }
    // Returns do not work with homogeneous coordinates, yet
};

/**
  * @private
  * 
  * Intersection of two circles using the stdform.
  * 
  * @param {Array} circ1 stdform of the first circle
  * @param {Array} circ2 stdform of the second circle
  * @param {number} i number of the returned intersection point. 
  *   i==0: use the positive square root, 
  *   i==1: use the negative square root.
  * @type JXG.Coords
  * @return Coordinates of the intersection point
  */
JXG.Algebra.prototype.meetCircleCircle = function(circ1,circ2,i) {
    var radicalAxis;
    if (circ1[4]<this.eps) { // Radius are zero, return center of circle, if on other circle
        if (this.distance(circ1.slice(1,3),circ2.slice(1,3))==circ2[4]) {
            return new JXG.Coords(JXG.COORDS_BY_USER, circ1.slice(1,3), this.board);
        } else {
            return new JXG.Coords(JXG.COORDS_BY_USER, [NaN,NaN], this.board);
        }
    }
    if (circ2[4]<this.eps) { // Radius are zero, return center of circle, if on other circle
        if (this.distance(circ2.slice(1,3),circ1.slice(1,3))==circ1[4]) {
            return new JXG.Coords(JXG.COORDS_BY_USER, circ2.slice(1,3), this.board);
        } else {
            return new JXG.Coords(JXG.COORDS_BY_USER, [NaN,NaN], this.board);
        }
    }
    radicalAxis = [circ2[3]*circ1[0]-circ1[3]*circ2[0],
                   circ2[3]*circ1[1]-circ1[3]*circ2[1],
                   circ2[3]*circ1[2]-circ1[3]*circ2[2],
                   0,1,Infinity, Infinity, Infinity];
    radicalAxis = this.normalize(radicalAxis);
    return this.meetLineCircle(radicalAxis,circ1,i);
    // Returns do not work with homogeneous coordinates, yet
};

/**
  * @private
  *
  * Normalize the stdform [c,b0,b1,a,k,r,q0,q1].
  * @param {Array} stdform to be normalized.
  * @type {Array}
  * @return The normalized stdform.
  */
JXG.Algebra.prototype.normalize = function(stdform) {
    var a2 = 2*stdform[3],
        r = stdform[4]/(a2),  // k/(2a)
        n, signr; 
    stdform[5] = r;
    stdform[6] = -stdform[1]/a2;
    stdform[7] = -stdform[2]/a2;
    if (r==Infinity || isNaN(r)) {
        n = Math.sqrt(stdform[1]*stdform[1]+stdform[2]*stdform[2]);
        stdform[0] /= n;
        stdform[1] /= n;
        stdform[2] /= n;
        stdform[3] = 0;
        stdform[4] = 1;
    } else if (Math.abs(r)>=1) {
        stdform[0] = (stdform[6]*stdform[6]+stdform[7]*stdform[7]-r*r)/(2*r);
        stdform[1] = -stdform[6]/r;
        stdform[2] = -stdform[7]/r;
        stdform[3] = 1/(2*r);
        stdform[4] = 1;
    } else {
        signr = (r<=0)?(-1):(1/*(r==0)?0:1*/);
        stdform[0] = signr*(stdform[6]*stdform[6]+stdform[7]*stdform[7]-r*r)*0.5;
        stdform[1] = -signr*stdform[6];
        stdform[2] = -signr*stdform[7];
        stdform[3] = signr/2;
        stdform[4] = signr*r;
    }
    return stdform;
};

/**
 * Compute an intersection of the curves c1 and c2
 * with a generalized Newton method.
 * We want to find values t1, t2 such that
 * c1(t1) = c2(t2), i.e.
 * (c1_x(t1)-c2_x(t2),c1_y(t1)-c2_y(t2)) = (0,0).
 * We set
 * (e,f) := (c1_x(t1)-c2_x(t2),c1_y(t1)-c2_y(t2))
 *
 * The Jacobian J is defined by
 * J = (a, b)
 *     (c, d)
 * where
 * a = c1_x'(t1)
 * b = -c2_x'(t2)
 * c = c1_y'(t1)
 * d = c2_y'(t2)
 *
 * The inverse J^(-1) of J is equal to
 *  (d, -b)/
 *  (-c, a) / (ad-bc)
 *
 * Then, (t1new, t2new) := (t1,t2) - J^(-1)*(e,f).
 * If the function meetCurveCurve possesses the properties
 * t1memo and t2memo then these are taken as start values
 * for the Newton algorithm.
 * After stopping of the Newton algorithm the values of t1 and t2 are stored in
 * t1memo and t2memo.
 * 
 * @param {JXG.Curve} c1: Curve, Line or Circle
 * @param {JXG.Curve} c2: Curve, Line or Circle
 * @param {float} t1ini: start value for t1
 * @param {float} t2ini: start value for t2
 * @type {JXG.Coords}
 * @return coordinate object for the intersection point
 **/
JXG.Algebra.prototype.meetCurveCurve = function(c1,c2,t1ini,t2ini) {
    var count = 0,
        t1, t2,
        a, b, c, d, disc,
        e, f, F, 
        D00, D01, 
        D10, D11;
        
    if (arguments.callee.t1memo) {
        t1 = arguments.callee.t1memo;
        t2 = arguments.callee.t2memo;
    } else {
        t1 = t1ini;
        t2 = t2ini;
    }
    if (t1>c1.maxX()) { t1 = c1.maxX(); }
    if (t1<c1.minX()) { t1 = c1.minX(); }
    if (t2>c2.maxX()) { t2 = c2.maxX(); }
    if (t2<c2.minX()) { t2 = c2.minX(); }
    e = c1.X(t1)-c2.X(t2);
    f = c1.Y(t1)-c2.Y(t2);
    F = e*e+f*f;
    
    D00 = c1.board.D(c1.X,c1);
    D01 = c2.board.D(c2.X,c2);
    D10 = c1.board.D(c1.Y,c1);
    D11 = c2.board.D(c2.Y,c2);
//$('debug').innerHTML = t1+' '+t2+'<br>\n';
    
    while (F>JXG.Math.eps && count<10) {
        a =  D00(t1);
        b = -D01(t2);
        c =  D10(t1);
        d = -D11(t2);
        disc = a*d-b*c;
        t1 -= (d*e-b*f)/disc;
        t2 -= (a*f-c*e)/disc;
        e = c1.X(t1)-c2.X(t2);
        f = c1.Y(t1)-c2.Y(t2);
        F = e*e+f*f;
        count++;
//$('debug').innerHTML += [a,b,c,d].join(':')+'['+disc+'], '+t1+' '+t2+ ' '+count+'<br>\n ';
    }

    arguments.callee.t1memo = t1;
    arguments.callee.t2memo = t2;
//    $('debug').innerHTML = arguments.callee.t1memo+' '+arguments.callee.t1memo+ ' '+count;
    //return (new JXG.Coords(JXG.COORDS_BY_USER, [2,2], this.board));
    if (Math.abs(t1)<Math.abs(t2)) {
        return (new JXG.Coords(JXG.COORDS_BY_USER, [c1.X(t1),c1.Y(t1)], this.board));
    } else {
        return (new JXG.Coords(JXG.COORDS_BY_USER, [c2.X(t2),c2.Y(t2)], this.board));
    }
};

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.

*/

/**
 * @fileoverview The object Intersection is defined in this file. Intersection
 * manages all properties and actiones required to cut circles and lines with
 * each other and draws the intersection points.
 * @author graphjs
 * @version 0.1
 */

/**
 * Constructs a new Intersection object.
 * @class This is the Intersection class. 
 * It manages all properties and actiones required to cut circles and lines with
 * each other and draws the intersection points.
 * @constructor
 * @param {String,Board} board The board the new point is drawn on.
 * @param {Array} coordinates An array with the affine user coordinates of the point.
 * @param {String} id Unique identifier for the point. If null or an empty string is given,
 *  an unique id will be generated by Board
 * @see JXG.Board#addPoint
 * @param {String} name Not necessarily unique name for the point. If null or an
 *  empty string is given, an unique name will be generated
 * @see JXG.Board#generateName
 * @param {bool} show False if the point is invisible, True otherwise
 */
JXG.Intersection = function(Board, Id, Intersect1, Intersect2, InterId1, InterId2, InterName1, InterName2) {
    this.constructor();
    /**
     * Reference to board where the intersected elements are drawn.
     * @type JXG.Board
     * @see JXG.Board
     */
    this.board = Board;
    
    /**
     * Unique identifier for the element. Equivalent to id-attribute of renderer element.
     * @type String
     */
    this.id = Id;
    this.name = this.id;

    /**
     * True when this object is visible, false otherwise.
     * @type bool
     */
    this.visProp = {};
    this.visProp['visible'] = true;
    this.show = true; // noch noetig? BV
    
    /**
     * True when the intersection points have real coordinates, false otherwise.
     * @type bool
     */    
    this.real = true;
    
    /** 
     * Stores all Intersection Objects which in this moment are not real and
     * hide this element.
     */
    this.notExistingParents = {};

    /**
     * Geometry element that is intersected with intersect2.
     * @type JXG.GeometryElement
     * @see #intersect2
     */
    this.intersect1 = JXG.getReference(this.board, Intersect1);

    /**
     * Geometry element that is intersected with intersect1.
     * @type JXG.GeometryElement
     * @see #intersect1
     */
    this.intersect2 = JXG.getReference(this.board, Intersect2);

    /**
     * Type of this object. For internal use only.
     * @private
     */
    this.type = JXG.OBJECT_TYPE_INTERSECTION;     

    /*
     * Only intersect existing geometry elements.
     */
    if( ((this.intersect1 == '') || (this.intersect1 == undefined)) && ((this.intersect2 == '') || (this.intersect2 == undefined))) {
        return;
    }

    /*
     * Do not intersect elements which aren't of type line, arrow, circle or arc.
     */
    if( ((this.intersect1.type == this.intersect2.type) && (this.intersect1.type == JXG.OBJECT_TYPE_LINE || this.intersect1.type == JXG.OBJECT_TYPE_ARROW)) 
         || ((this.intersect1.type == JXG.OBJECT_TYPE_LINE) && (this.intersect2.type == JXG.OBJECT_TYPE_ARROW))
         || ((this.intersect2.type == JXG.OBJECT_TYPE_LINE) && (this.intersect1.type == JXG.OBJECT_TYPE_ARROW)) ) {
        /* Intersect two elements of type line or arrow */
        
        var coords = this.board.algebra.intersectLineLine(this.intersect1, this.intersect2).usrCoords.slice(1);

        /* Create intersection point */
        this.p = new JXG.Point(this.board, coords, InterId1, InterName1, true);
        /* A point constructed by an intersection can't be moved, so it is fixed */
        this.p.fixed = true;
        this.addChild(this.p);
        this.real = true;

        /* 
         * Because the update function depends on the types of the intersected elements
         * the update method has to be defined dynamically in dependence of the intersected
         * elements.
         */
        this.update = function () {
            /* Calculate the coordinates of the intersection point in dependance of the intersected elements */
            if (this.needsUpdate) {
                this.p.coords = this.board.algebra.intersectLineLine(this.intersect1, this.intersect2);
                /* Update the point */
                //this.p.update();
                this.needsUpdate = false;
            }
        };
        
        /*
         * Hides the element, generated dynamically.
         */
        this.hideElement = function() {
            this.visProp['visible'] = false;
            this.p.hideElement();
        };
        
        /*
         * Shows the element, generated dynamically.
         */
        this.showElement = function() {
            this.visProp['visible'] = true;
            this.p.showElement();
        };
        
        /*
         * Hides the element and his children. This is called from parents which became invisible or unreal
         * and so this element isn't real anymore. The not existing parent is stored in the notExistingParents
         * array.
         */
        this.hideChild = function(id) {    
            this.notExistingParents[id] = this.board.objects[id];

            for(var el in this.descendants) {    
                if(this.descendants[el].visProp['visible'] && this.descendants[el].type != JXG.OBJECT_TYPE_INTERSECTION) {
                    if(this.descendants[el].type != JXG.OBJECT_TYPE_TEXT) {
                        this.descendants[el].hideElement();
                        this.descendants[el].visProp['visible'] = true;
                    }
                    else {
                        if(!this.descendants[el].isLabel) {
                            this.descendants[el].hideElement();
                            this.descendants[el].visProp['visible'] = true;
                        }
                    }
                }      
                this.descendants[el].notExistingParents[id] = this.board.objects[id];
            }            
        };

        /*
         * Shows the element and his children. This is called from parents which became visible or real
         * and so this element is now real. The formerly not existing parent is deleted from the
         * notExistingParents array.
         */
        this.showChild = function(id) {        
            for(var el in this.board.objects) {        
                delete(this.board.objects[el].notExistingParents[id]);
                if(this.board.objects[el].visProp['visible'] && JXG.keys(this.board.objects[el].notExistingParents).length == 0) {
                    if(this.board.objects[el].type != JXG.OBJECT_TYPE_INTERSECTION) {
                        this.board.objects[el].showElement();
                    }
                }
            }
        };
    }
    else if( ((Intersect1.type == Intersect2.type) && (Intersect1.type == JXG.OBJECT_TYPE_CIRCLE || Intersect1.type == JXG.OBJECT_TYPE_ARC)) ||
              (Intersect1.type == JXG.OBJECT_TYPE_CIRCLE && Intersect2.type == JXG.OBJECT_TYPE_ARC) ||
              (Intersect2.type == JXG.OBJECT_TYPE_CIRCLE && Intersect1.type == JXG.OBJECT_TYPE_ARC) ) { // Circle <-> Circle, Arc <-> Arc, Arc <-> Circle,
        this.p1 = new JXG.Point(this.board, [0, 0], InterId1, InterName1, false);
        this.p1.fixed = true;
        this.p1.label.content.visProp['visible'] = true;
        this.p2 = new JXG.Point(this.board, [0, 0], InterId2, InterName2, false);
        this.p2.fixed = true;
        this.p2.label.content.visProp['visible'] = true;
        this.addChild(this.p1);
        this.addChild(this.p2);

        var coordinates = this.board.algebra.intersectCircleCircle(this.intersect1, this.intersect2);
        if(coordinates[0] == 1) {
            this.p1.coords = coordinates[1];
            this.p1.showElement();
            this.p1.updateRenderer();

            this.p2.coords = coordinates[2];
            this.p2.showElement();
            this.p2.updateRenderer();
            
            this.real = true;
        }
        else {
            this.real = false;
        }

        this.update = function () {    
            if (!this.needsUpdate) { return; }
            var coordinates = this.board.algebra.intersectCircleCircle(this.intersect1, this.intersect2);
            var p1show = this.p1.visProp['visible'];
            var p2show = this.p2.visProp['visible'];         
            if(coordinates[0] == 0) {  
                if(this.real) {
                    this.hideChild(this.id);
                    this.p1.visProp['visible'] = p1show;
                    this.p2.visProp['visible'] = p2show;
                    this.real = false;
                }
            } else {
                this.p1.coords = coordinates[1];     
                this.p2.coords = coordinates[2];
                if(!this.real) {
                    this.showChild(this.id); 
                    this.real = true;
                }
            }
            this.needsUpdate = false;
        };
        
        this.hideElement = function() {
            this.visProp['visible'] = false;
            this.p1.hideElement();
            this.p2.hideElement();
        };
        
        this.showElement = function() {
            this.visProp['visible'] = true;
            this.p1.showElement();
            this.p2.showElement();
        };

        this.hideChild = function(id) {
            this.notExistingParents[id] = this.board.objects[id];

            for(var el in this.descendants) {    
                if(this.descendants[el].visProp['visible'] && this.descendants[el].type != JXG.OBJECT_TYPE_INTERSECTION) {
                    if(this.descendants[el].type != JXG.OBJECT_TYPE_TEXT) {
                        this.descendants[el].hideElement();
                        this.descendants[el].visProp['visible'] = true;
                    }
                    else {
                        if(!this.descendants[el].isLabel) {
                            this.descendants[el].hideElement();
                            this.descendants[el].visProp['visible'] = true;
                        }
                    }
                }      
                this.descendants[el].notExistingParents[id] = this.board.objects[id];
            }                
        };

        this.showChild = function(id) {            
            for(el in this.board.objects) {        
                delete(this.board.objects[el].notExistingParents[id]);
                if(this.board.objects[el].visProp['visible'] && JXG.keys(this.board.objects[el].notExistingParents).length == 0) {
                    if(this.board.objects[el].type != JXG.OBJECT_TYPE_INTERSECTION) {
                        this.board.objects[el].showElement();
                    }
                }        
            }            
        };
    }
    else { // Circle <-> Line, Arc <-> Line, Circle <-> Arrow, Arc <-> Arrow
        this.p1 = new JXG.Point(this.board, [0, 0], InterId1, InterName1, false);
        this.p1.fixed = true;
        this.p1.label.content.visProp['visible'] = true;        
        this.p2 = new JXG.Point(this.board, [0, 0], InterId2, InterName2, false);
        this.p2.fixed = true;
        this.p2.label.content.visProp['visible'] = true;
        this.addChild(this.p1);
        this.addChild(this.p2);        
        
        if(this.intersect1.type == JXG.OBJECT_TYPE_LINE || this.intersect1.type == JXG.OBJECT_TYPE_ARROW) {
            var swap = this.intersect1;
            this.intersect1 = this.intersect2;
            this.intersect2 = swap;
        }
        
        var coordinates = this.board.algebra.intersectCircleLine(this.intersect1, this.intersect2);
        if(coordinates[0] == 1) { // not really implemented
            this.p1.coords = coordinates[1];
            this.p1.showElement();
            this.p1.update();    
        } 
        else if(coordinates[0] == 2) {
            this.p1.coords = coordinates[1];
            this.p1.showElement();        

            this.p2.coords = coordinates[2];
            this.p2.showElement();            

            //this.p1.update();
            this.p1.updateRenderer();
            //this.p2.update(); 
            this.p2.updateRenderer();    
            
            this.real = true;
        }
        else {
            this.real = false;
        }

        this.update = function () {
            if (!this.needsUpdate) { return; }
            var coordinates = this.board.algebra.intersectCircleLine(this.intersect1, this.intersect2);
            var show1 = this.p1.visProp['visible'];
            var show2 = this.p2.visProp['visible'];
            
            if(coordinates[0] == 0) {
                if(this.real) {
                    this.hideChild(this.id);
                    this.p1.visProp['visible'] = show1; 
                    this.p2.visProp['visible'] = show2;
                    this.real = false;
                }
            } else if(coordinates[0] == 2) {
                this.p1.coords = coordinates[1];   
                this.p2.coords = coordinates[2];
                if(!this.real) {
                    this.showChild(this.id);  
                    this.real = true;
                }
            }
            this.needsUpdate = false;
        };
        
        this.hideElement = function() {
            this.visProp['visible'] = false;
            this.p1.hideElement();
            this.p2.hideElement();
        };
        
        this.showElement = function() {
            this.visProp['visible'] = true;
            this.p1.showElement();
            this.p2.showElement();
        };

        this.hideChild = function(id) {
            this.notExistingParents[id] = this.board.objects[id];

            for(var el in this.descendants) {    
                if(this.descendants[el].visProp['visible'] && this.descendants[el].type != JXG.OBJECT_TYPE_INTERSECTION) {
                    if(this.descendants[el].type != JXG.OBJECT_TYPE_TEXT) {
                        this.descendants[el].hideElement();
                        this.descendants[el].visProp['visible'] = true;
                    }
                    else {
                        if(!this.descendants[el].isLabel) {
                            this.descendants[el].hideElement();
                            this.descendants[el].visProp['visible'] = true;
                        }
                    }
                }      
                this.descendants[el].notExistingParents[id] = this.board.objects[id];
            }                
        };

        this.showChild = function(id) {
            for(el in this.board.objects) {            
                delete(this.board.objects[el].notExistingParents[id]);
                if(this.board.objects[el].visProp['visible'] && JXG.keys(this.board.objects[el].notExistingParents).length == 0) {
                    if(this.board.objects[el].type != JXG.OBJECT_TYPE_INTERSECTION) {
                        this.board.objects[el].showElement();
                    }
                }            
            }            
        };
    }

    this.id = this.board.addIntersection(this);
};
JXG.Intersection.prototype = new JXG.GeometryElement();
   
/**
 * Calls the renderer to update the drawing. This method is defined dynamically
 * as it highly depends on the types of the intersected elements.
 */
JXG.Intersection.prototype.update = function() {
    return;
};

/**
 * Checks whether (x,y) is near the point.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @return {bool} Always returns false
 */
JXG.Intersection.prototype.hasPoint = function(x, y) {
    return false;
};

/**
 * Hides the element and his children. This is called from parents which became invisible or unreal
 * and so this element isn't real anymore. The not existing parent is stored in the notExistingParents
 * array.
 * @param {String} id The identifier of the element causing this element to be hidden.
 */
JXG.Intersection.prototype.hideChild = function(id) {
};

/**
 * Shows the element and his children. This is called from parents which became visible or real
 * and so this element is now real. The formerly not existing parent is deleted from the
 * notExistingParents array.
 * @param {String} id The identifier of the element causing this element to be shown.
 */
JXG.Intersection.prototype.showChild = function(id) {
};

/**
 * Remove intersection points from drawing.
 */
JXG.Intersection.prototype.remove = function() {
    if(this.p != undefined)
        this.board.removeObject(this.p);
    if(this.p1 != undefined)
        this.board.removeObject(this.p1);
    if(this.p2 != undefined)
        this.board.removeObject(this.p2);
        
    return;
};

/**
 * Dummy method 
 */
JXG.Intersection.prototype.updateRenderer = function() {
};

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview This file contains our composition elements, i.e. these elements are mostly put together
 * from one or more {@link JXG.GeometryElement} but with a special meaning. E.g. the midpoint element is contained here
 * and this is just a {@link JXG.Point} with coordinates dependent from two other points. Currently in this file the
 * following compositions can be found: <ul>
 *   <li>{@link Arrowparallel} (currently private)</li>
 *   <li>{@link Bisector}</li>
 *   <li>{@link Circumcircle}</li>
 *   <li>{@link Circumcirclemidpoint}</li>
 *   <li>{@link Integral}</li>
 *   <li>{@link Midpoint}</li>
 *   <li>{@link Mirrorpoint}</li>
 *   <li>{@link Normal}</li>
 *   <li>{@link Parallel}</li>
 *   <li>{@link Perpendicular}</li>
 *   <li>{@link Perpendicularpoint}</li>
 *   <li>{@link Reflection}</li></ul>
 */

/**
 * @class This is used to construct a perpendicular point.
 * @pseudo
 * @description A perpendicular point is given by a point and a line. It is determined by projecting the given point
 * orthogonal onto the given line.
 * @constructor
 * @name Perpendicularpoint
 * @type JXG.Point
 * @augments JXG.Point
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Line_JXG.Point} p,l The constructed point is the orthogonal projection of p onto l.
 * @example
 * var p1 = board.create('point', [0.0, 4.0]);
 * var p2 = board.create('point', [6.0, 1.0]);
 * var l1 = board.create('line', [p1, p2]);
 * var p3 = board.create('point', [3.0, 3.0]);
 *
 * var pp1 = board.create('perpendicularpoint', [p3, l1]);
 * </pre><div id="ded148c9-3536-44c0-ab81-1bb8fa48f3f4" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var ppex1_board = JXG.JSXGraph.initBoard('ded148c9-3536-44c0-ab81-1bb8fa48f3f4', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var ppex1_p1 = ppex1_board.create('point', [0.0, 4.0]);
 *   var ppex1_p2 = ppex1_board.create('point', [6.0, 1.0]);
 *   var ppex1_l1 = ppex1_board.create('line', [ppex1_p1, ppex1_p2]);
 *   var ppex1_p3 = ppex1_board.create('point', [3.0, 3.0]);
 *   var ppex1_pp1 = ppex1_board.create('perpendicularpoint', [ppex1_p3, ppex1_l1]);
 * </script><pre>
 */
JXG.createPerpendicularPoint = function(board, parentArr, atts) {
    var l, p, t;

    if(JXG.isPoint(parentArr[0]) && parentArr[1].type == JXG.OBJECT_TYPE_LINE) {
        p = parentArr[0];
        l = parentArr[1];
    }
    else if(JXG.isPoint(parentArr[1]) && parentArr[0].type == JXG.OBJECT_TYPE_LINE) {
        p = parentArr[1];
        l = parentArr[0];
    }
    else {
        throw new Error("JSXGraph: Can't create perpendicular point with parent types '" + (typeof parentArr[0]) + "' and '" + (typeof parentArr[1]) + "'.");
    }

    // no need to call create, the properties will be set through the create('perpendicular') call
    t = JXG.createPoint(board, [function () { return board.algebra.perpendicular(l, p)[0]; }], {fixed: true, name: atts['name'], id: atts['id']});
    p.addChild(t); // notwendig, um auch den Punkt upzudaten
    l.addChild(t);

    t.update();

    t.generatePolynomial = function() {
        /*
         *  Perpendicular takes point P and line L and creates point T and line M:
         *
         *                          | M
         *                          |
         *                          x P (p1,p2)
         *                          |
         *                          |
         *  L                       |
         *  ----------x-------------x------------------------x--------
         *            A (a1,a2)     |T (t1,t2)               B (b1,b2)
         *                          |
         *                          |
         *
         * So we have two conditions:
         *
         *   (a)  AT  || TB          (collinearity condition)
         *   (b)  PT _|_ AB          (orthogonality condition)
         *
         *      a2-t2       t2-b2
         *     -------  =  -------           (1)
         *      a1-t1       t1-b1
         *
         *      p2-t2         a1-b1
         *     -------  =  - -------         (2)
         *      p1-t1         a2-b2
         *
         * Multiplying (1) and (2) with denominators and simplifying gives
         *
         *    a2t1 - a2b1 + t2b1 - a1t2 + a1b2 - t1b2 = 0                  (1')
         *
         *    p2a2 - p2b2 - t2a2 + t2b2 + p1a1 - p1b1 - t1a1 + t1b1 = 0    (2')
         *
         */

        var a1 = l.point1.symbolic.x;
        var a2 = l.point1.symbolic.y;
        var b1 = l.point2.symbolic.x;
        var b2 = l.point2.symbolic.y;
        var p1 = p.symbolic.x;
        var p2 = p.symbolic.y;
        var t1 = t.symbolic.x;
        var t2 = t.symbolic.y;

        var poly1 = '('+a2+')*('+t1+')-('+a2+')*('+b1+')+('+t2+')*('+b1+')-('+a1+')*('+t2+')+('+a1+')*('+b2+')-('+t1+')*('+b2+')';
        var poly2 = '('+p2+')*('+a2+')-('+p2+')*('+b2+')-('+t2+')*('+a2+')+('+t2+')*('+b2+')+('+p1+')*('+a1+')-('+p1+')*('+b1+')-('+t1+')*('+a1+')+('+t1+')*('+b1+')';

        return [poly1, poly2];
    };

    return t;
};


/**
 * @class This element is used to provide a constructor for a perpendicular.
 * @pseudo
 * @description  A perpendicular is a composition of two elements: a line and a point. The line is orthogonal
 * to a given line and contains a given point and meets the given line in the perpendicular point.
 * @name Perpendicular
 * @constructor
 * @type array
 * @return An array containing two elements: A {@link JXG.Line} object in the first component and a
 * {@link JXG.Point} element in the second component. The line is orthogonal to the given line and meets it
 * in the returned point.
 * @throws {Exception} If the elements cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Line_JXG.Point} l,p The perpendicular line will be orthogonal to l and
 * will contain p. The perpendicular point is the intersection point of the two lines.
 * @example
 * // Create a perpendicular
 * var p1 = board.create('point', [0.0, 2.0]);
 * var p2 = board.create('point', [2.0, 1.0]);
 * var l1 = board.create('line', [p1, p2]);
 *
 * var p3 = board.create('point', [3.0, 3.0]);
 * var perp1 = board.create('perpendicular', [l1, p3]);
 * </pre><div id="037a6eb2-781d-4b71-b286-763619a63f22" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var pex1_board = JXG.JSXGraph.initBoard('037a6eb2-781d-4b71-b286-763619a63f22', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var pex1_p1 = pex1_board.create('point', [0.0, 2.0]);
 *   var pex1_p2 = pex1_board.create('point', [2.0, 1.0]);
 *   var pex1_l1 = pex1_board.create('line', [pex1_p1, pex1_p2]);
 *   var pex1_p3 = pex1_board.create('point', [3.0, 3.0]);
 *   var pex1_perp1 = pex1_board.create('perpendicular', [pex1_l1, pex1_p3]);
 * </script><pre>
 */
JXG.createPerpendicular = function(board, parentArr, atts) {
    var p, l, pd, t, ret;

    parentArr[0] = JXG.getReference(board, parentArr[0]);
    parentArr[1] = JXG.getReference(board, parentArr[1]);

    if(JXG.isPoint(parentArr[0]) && parentArr[1].elementClass == JXG.OBJECT_CLASS_LINE) {
        l = parentArr[1];
        p = parentArr[0];
    }
    else if(JXG.isPoint(parentArr[1]) && parentArr[0].elementClass == JXG.OBJECT_CLASS_LINE) {
        l = parentArr[0];
        p = parentArr[1];
    }
    else {
        throw new Error("JSXGraph: Can't create perpendicular with parent types '" + (typeof parentArr[0]) + "' and '" + (typeof parentArr[1]) + "'.");
    }

    if(!JXG.isArray(atts['id'])) {
        atts['id'] = ['',''];
    }
    if(!JXG.isArray(atts['name'])) {
        atts['name'] = ['',''];
    }

    // no need to call create, the properties will be set through the create('perpendicular') call
    t = JXG.createPerpendicularPoint(board, [l, p], {fixed: true, name: atts['name'][1], id: atts['id'][1], visible: false});
    pd = JXG.createSegment(board, [function () { return (board.algebra.perpendicular(l, p)[1] ? [t, p] : [p, t]); }], {name: atts['name'][0], id: atts['id'][0]});

    ret = [pd, t];
    ret.line = pd;
    ret.point = t;
    ret.multipleElements = true;

    return ret;
};

/**
 * @class The midpoint element constructs a point in the middle of two given points.
 * @pseudo
 * @description A midpoint is given by two points. It is collinear to the given points and the distance
 * is the same to each of the given points, i.e. it is in the middle of the given points.
 * @constructor
 * @name Midpoint
 * @type JXG.Point
 * @augments JXG.Point
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point_JXG.Point} p1,p2 The constructed point will be in the middle of p1 and p2.
 * @param {JXG.Line} l The midpoint will be in the middle of {@link JXG.Line#point1} and {@link JXG.Line#point2} of
 * the given line l.
 * @example
 * // Create base elements: 2 points and 1 line
 * var p1 = board.create('point', [0.0, 2.0]);
 * var p2 = board.create('point', [2.0, 1.0]);
 * var l1 = board.create('segment', [[0.0, 3.0], [3.0, 3.0]]);
 *
 * var mp1 = board.create('midpoint', [p1, p2]);
 * var mp2 = board.create('midpoint', [l1]);
 * </pre><div id="7927ef86-24ae-40cc-afb0-91ff61dd0de7" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var mpex1_board = JXG.JSXGraph.initBoard('7927ef86-24ae-40cc-afb0-91ff61dd0de7', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var mpex1_p1 = mpex1_board.create('point', [0.0, 2.0]);
 *   var mpex1_p2 = mpex1_board.create('point', [2.0, 1.0]);
 *   var mpex1_l1 = mpex1_board.create('segment', [[0.0, 3.0], [3.0, 3.0]]);
 *   var mpex1_mp1 = mpex1_board.create('midpoint', [mpex1_p1, mpex1_p2]);
 *   var mpex1_mp2 = mpex1_board.create('midpoint', [mpex1_l1]);
 * </script><pre>
 */
JXG.createMidpoint = function(board, parentArr, atts) {
    var a, b, t;
    if(parentArr.length == 2 && JXG.isPoint(parentArr[0]) && JXG.isPoint(parentArr[1])) {
        a = parentArr[0];
        b = parentArr[1];
    }
    else if(parentArr.length == 1 && parentArr[0].elementClass == JXG.OBJECT_CLASS_LINE) {
        a = parentArr[0].point1;
        b = parentArr[0].point2;
    }
    else {
        throw new Error("JSXGraph: Can't create midpoint.");
    }

    if(atts) {
    	atts['fixed'] = true;
    } else {
    	atts = {fixed: true};
    }

    t = board.create('point', [function () { return (a.coords.usrCoords[1] + b.coords.usrCoords[1])/2.; },
                                      function () { return (a.coords.usrCoords[2] + b.coords.usrCoords[2])/2.; }], atts);
    a.addChild(t);
    b.addChild(t);

    t.update();

    t.generatePolynomial = function() {
        /*
         *  Midpoint takes two point A and B or line L (with points P and Q) and creates point T:
         *
         *  L (not necessarily)
         *  ----------x------------------x------------------x--------
         *            A (a1,a2)          T (t1,t2)          B (b1,b2)
         *
         * So we have two conditions:
         *
         *   (a)   AT  ||  TB           (collinearity condition)
         *   (b)  [AT] == [TB]          (equidistant condition)
         *
         *      a2-t2       t2-b2
         *     -------  =  -------                                         (1)
         *      a1-t1       t1-b1
         *
         *     (a1 - t1)^2 + (a2 - t2)^2 = (b1 - t1)^2 + (b2 - t2)^2       (2)
         *
         *
         * Multiplying (1) with denominators and simplifying (1) and (2) gives
         *
         *    a2t1 - a2b1 + t2b1 - a1t2 + a1b2 - t1b2 = 0                      (1')
         *
         *    a1^2 - 2a1t1 + a2^2 - 2a2t2 - b1^2 + 2b1t1 - b2^2 + 2b2t2 = 0    (2')
         *
         */

        var a1 = a.symbolic.x;
        var a2 = a.symbolic.y;
        var b1 = b.symbolic.x;
        var b2 = b.symbolic.y;
        var t1 = t.symbolic.x;
        var t2 = t.symbolic.y;

        var poly1 = '('+a2+')*('+t1+')-('+a2+')*('+b1+')+('+t2+')*('+b1+')-('+a1+')*('+t2+')+('+a1+')*('+b2+')-('+t1+')*('+b2+')';
        var poly2 = '('+a1+')^2 - 2*('+a1+')*('+t1+')+('+a2+')^2-2*('+a2+')*('+t2+')-('+b1+')^2+2*('+b1+')*('+t1+')-('+b2+')^2+2*('+b2+')*('+t2+')';

        return [poly1, poly2];
    };

    // 2009/10/11, mg:
    // * this is just a test: we're forming three closures in here: One with generatePolynomial and two when we create the
    //   point by using function objects as positions. So a reference to the current Activation object is held in the scope of
    //   those three functions. But we don't need a reference to the atts object in any of them, so we're deleting the
    //   reference to it. Hence, the garbage collector can remove it from memory, if there's no reference to it left.
    // * first test successful: attributes will be set properly. it remains to test, if this helps us to reduce memory usage.
    //   for the test we'll query a JXG attribute "nullAtts" which has to be set to true to null the atts parameter. Then
    //   50,000 midpoints will be created in an example file with resp. without nulling the atts parameter and after that chrome
    //   will be used to compare the memory usage.
    if(JXG.nullAtts)
        atts = null;

    return t;
};

/**
 * @class This element is used to construct a parallel point.
 * @pseudo
 * @description A parallel point is given by three points. Taking the euclidean vector from the first to the
 * second point, the parallel point is determined by adding that vector to the third point.
 * The line determined by the first two points is parallel to the line determined by the third point and the constructed point.
 * @constructor
 * @name Parallelpoint
 * @type JXG.Point
 * @augments JXG.Point
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point_JXG.Point_JXG.Point} p1,p2,p3 Taking the euclidean vector <tt>v=p2-p1</tt> the parallel point is determined by
 * <tt>p4 = p3+v</tt>
 * @param {JXG.Line_JXG.Point} l,p The resulting point will together with p specify a line which is parallel to l.
 * @example
 * var p1 = board.create('point', [0.0, 2.0]);
 * var p2 = board.create('point', [2.0, 1.0]);
 * var p3 = board.create('point', [3.0, 3.0]);
 *
 * var pp1 = board.create('parallelpoint', [p1, p2, p3]);
 * </pre><div id="488c4be9-274f-40f0-a469-c5f70abe1f0e" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var ppex1_board = JXG.JSXGraph.initBoard('488c4be9-274f-40f0-a469-c5f70abe1f0e', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var ppex1_p1 = ppex1_board.create('point', [0.0, 2.0]);
 *   var ppex1_p2 = ppex1_board.create('point', [2.0, 1.0]);
 *   var ppex1_p3 = ppex1_board.create('point', [3.0, 3.0]);
 *   var ppex1_pp1 = ppex1_board.create('parallelpoint', [ppex1_p1, ppex1_p2, ppex1_p3]);
 * </script><pre>
 */
JXG.createParallelPoint = function(board, parentArr, atts) {
    var a, b, c, p;

    if(parentArr.length == 3 && parentArr[0].elementClass == JXG.OBJECT_CLASS_POINT && parentArr[1].elementClass == JXG.OBJECT_CLASS_POINT && parentArr[2].elementClass == JXG.OBJECT_CLASS_POINT) {
        a = parentArr[0];
        b = parentArr[1];
        c = parentArr[2];
    } else if (parentArr[0].elementClass == JXG.OBJECT_CLASS_POINT && parentArr[1].elementClass == JXG.OBJECT_CLASS_LINE) {
        c = parentArr[0];
        a = parentArr[1].point1;
        b = parentArr[1].point2;
    } else if (parentArr[1].elementClass == JXG.OBJECT_CLASS_POINT && parentArr[0].elementClass == JXG.OBJECT_CLASS_LINE) {
        c = parentArr[1];
        a = parentArr[0].point1;
        b = parentArr[0].point2;
    }
    else {
        throw new Error("JSXGraph: Can't create parallel point with parent types '" + (typeof parentArr[0]) + "', '" + (typeof parentArr[1]) + "' and '" + (typeof parentArr[2]) + "'.");
    }

    p = board.create('point', [function () { return c.coords.usrCoords[1] + b.coords.usrCoords[1] - a.coords.usrCoords[1]; }, function () { return c.coords.usrCoords[2] + b.coords.usrCoords[2] - a.coords.usrCoords[2]; }], atts);
//    p1.addChild(p); // required for algorithms requiring dependencies between elements
//    p2.addChild(p);
    c.addChild(p);

    // required to set the coordinates because functions are considered as constraints. hence, the coordinates get set first after an update.
    // can be removed if the above issue is resolved.
    p.update();

    p.generatePolynomial = function() {
        /*
         *  Parallelpoint takes three points A, B and C or line L (with points B and C) and creates point T:
         *
         *
         *                     C (c1,c2)                             T (t1,t2)
         *                      x                                     x
         *                     /                                     /
         *                    /                                     /
         *                   /                                     /
         *                  /                                     /
         *                 /                                     /
         *                /                                     /
         *               /                                     /
         *              /                                     /
         *  L (opt)    /                                     /
         *  ----------x-------------------------------------x--------
         *            A (a1,a2)                             B (b1,b2)
         *
         * So we have two conditions:
         *
         *   (a)   CT  ||  AB           (collinearity condition I)
         *   (b)   BT  ||  AC           (collinearity condition II)
         *
         * The corresponding equations are
         *
         *    (b2 - a2)(t1 - c1) - (t2 - c2)(b1 - a1) = 0         (1)
         *    (t2 - b2)(a1 - c1) - (t1 - b1)(a2 - c2) = 0         (2)
         *
         * Simplifying (1) and (2) gives
         *
         *    b2t1 - b2c1 - a2t1 + a2c1 - t2b1 + t2a1 + c2b1 - c2a1 = 0      (1')
         *    t2a1 - t2c1 - b2a1 + b2c1 - t1a2 + t1c2 + b1a2 - b1c2 = 0      (2')
         *
         */

        var a1 = a.symbolic.x;
        var a2 = a.symbolic.y;
        var b1 = b.symbolic.x;
        var b2 = b.symbolic.y;
        var c1 = c.symbolic.x;
        var c2 = c.symbolic.y;
        var t1 = p.symbolic.x;
        var t2 = p.symbolic.y;

        var poly1 =  '('+b2+')*('+t1+')-('+b2+')*('+c1+')-('+a2+')*('+t1+')+('+a2+')*('+c1+')-('+t2+')*('+b1+')+('+t2+')*('+a1+')+('+c2+')*('+b1+')-('+c2+')*('+a1+')';
        var poly2 =  '('+t2+')*('+a1+')-('+t2+')*('+c1+')-('+b2+')*('+a1+')+('+b2+')*('+c1+')-('+t1+')*('+a2+')+('+t1+')*('+c2+')+('+b1+')*('+a2+')-('+b1+')*('+c2+')';

        return [poly1, poly2];
    };

    return p;
};

/**
 * @class Constructor for a parallel line.
 * @pseudo
 * @description A parallel is a line through a given point with the same slope as a given line.
 * @constructor
 * @name Parallel
 * @type JXG.Line
 * @augments JXG.Line
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Line_JXG.Point} l,p The constructed line contains p and has the same slope as l.
 * @example
 * // Create a parallel
 * var p1 = board.create('point', [0.0, 2.0]);
 * var p2 = board.create('point', [2.0, 1.0]);
 * var l1 = board.create('line', [p1, p2]);
 *
 * var p3 = board.create('point', [3.0, 3.0]);
 * var pl1 = board.create('parallel', [l1, p3]);
 * </pre><div id="24e54f9e-5c4e-4afb-9228-0ef27a59d627" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var plex1_board = JXG.JSXGraph.initBoard('24e54f9e-5c4e-4afb-9228-0ef27a59d627', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var plex1_p1 = plex1_board.create('point', [0.0, 2.0]);
 *   var plex1_p2 = plex1_board.create('point', [2.0, 1.0]);
 *   var plex1_l1 = plex1_board.create('line', [plex1_p1, plex1_p2]);
 *   var plex1_p3 = plex1_board.create('point', [3.0, 3.0]);
 *   var plex1_pl1 = plex1_board.create('parallel', [plex1_l1, plex1_p3]);
 * </script><pre>
 */
JXG.createParallel = function(board, parents, atts) {
    var p, pp, pl, cAtts;

    /* parallel point polynomials are done in createParallelPoint */

    cAtts = {name: null, id: null, fixed: true, visible: false};
    if(JXG.isArray(atts['name']) && atts['name'].length == 2) {
        cAtts['name'] = atts['name'][1];
        atts['name'] = atts['name'][0];
    } else
        cAtts['name'] = atts['name'] + 'p2';
    if(JXG.isArray(atts['id']) && atts['id'].length == 2) {
        cAtts['id'] = atts['id'][1];
        atts['id'] = atts['id'][0];
    } else
        cAtts['id'] = atts['id'] + 'p2';

    if(atts) {
        cAtts = JXG.cloneAndCopy(atts, cAtts);
    }

    try {
        pp = JXG.createParallelPoint(board, parents, cAtts);
    } catch (e) {
        throw new Error("JSXGraph: Can't create parallel with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
    }

    p = null;
    if(parents.length == 3)
        p = parents[2];
    else if (parents[0].elementClass == JXG.OBJECT_CLASS_POINT)
        p = parents[0];
    else if (parents[1].elementClass == JXG.OBJECT_CLASS_POINT)
        p = parents[1];

    pl = board.create('line', [p, pp], atts);

    return pl;
};

/**
 * TODO is this really required? it is the same as 'parallel', except that it doesn't touch the first/lastarrow properties and it returns
 * the parallel point. for now it is set to private. please review the docs-comment before making it public. especially the example section
 * isn't done by now. --michael
 * @private
 * @class Constructs two elements: an arrow and a point.
 * @pseudo
 * @description An arrow parallel is an arrow through a given point with the same slope as another given arrow.
 * @constructor
 * @name Arrowparallel
 * @type JXG.Line
 * @augments JXG.Line
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {Arrow_JXG.Point} a,p The constructed arrow contains p and has the same slope as a.
 * @example
 * // Create a parallel
 * var p1 = board.create('point', [0.0, 2.0]);
 * var p2 = board.create('point', [2.0, 1.0]);
 * var l1 = board.create('line', [p1, p2]);
 *
 * var p3 = board.create('point', [3.0, 3.0]);
 * var pl1 = board.create('parallel', [l1, p3]);
 * </pre><div id="qwe" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var plex1_board = JXG.JSXGraph.initBoard('asd', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var plex1_p1 = plex1_board.create('point', [0.0, 2.0]);
 *   var plex1_p2 = plex1_board.create('point', [2.0, 1.0]);
 *   var plex1_l1 = plex1_board.create('line', [plex1_p1, plex1_p2]);
 *   var plex1_p3 = plex1_board.create('point', [3.0, 3.0]);
 *   var plex1_pl1 = plex1_board.create('parallel', [plex1_l1, plex1_p3]);
 * </script><pre>
 */
JXG.createArrowParallel = function(board, parents, atts) {
    var l, cAtts;

    /* parallel arrow point polynomials are done in createParallelPoint */
    try {
        // we don't have to get onto that whole create stack here
        // because that'll be run for the line l right after leaving that function.
        l = JXG.createParallel(board, parents, atts);
    } catch (e) {
        throw new Error("JSXGraph: Can't create arrowparallel with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
    }

    // Select the default behavior. If the user wants something else he would set it in atts.
    // That gets parsed and set right after this function.
    l.setStraight(false, false);
    l.setArrow(false,true);
    return l;
};

/**
 * @class Constructs a normal.
 * @pseudo
 * @description A normal is a line through a given point on a element of type line, circle, curve, or turtle and orthogonal to that object.
 * @constructor
 * @name Normal
 * @type JXG.Line
 * @augments JXG.Line
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Line,JXG.Circle,JXG.Curve,JXG.Turtle_JXG.Point} o,p The constructed line contains p which lies on the object and is orthogonal
 * to the tangent to the object in the given point.
 * @param {Glider} p Works like above, however the object is given by {@link Glider#slideObject}.
 * @example
 * // Create a normal to a circle.
 * var p1 = board.create('point', [2.0, 2.0]);
 * var p2 = board.create('point', [3.0, 2.0]);
 * var c1 = board.create('circle', [p1, p2]);
 *
 * var norm1 = board.create('normal', [c1, p2]);
 * </pre><div id="4154753d-3d29-40fb-a860-0b08aa4f3743" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var nlex1_board = JXG.JSXGraph.initBoard('4154753d-3d29-40fb-a860-0b08aa4f3743', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var nlex1_p1 = nlex1_board.create('point', [2.0, 2.0]);
 *   var nlex1_p2 = nlex1_board.create('point', [3.0, 2.0]);
 *   var nlex1_c1 = nlex1_board.create('circle', [nlex1_p1, nlex1_p2]);
 *
 *   // var nlex1_p3 = nlex1_board.create('point', [1.0, 2.0]);
 *   var nlex1_norm1 = nlex1_board.create('normal', [nlex1_c1, nlex1_p2]);
 * </script><pre>
 */
JXG.createNormal = function(board, parents, attributes) {
    /* TODO normal polynomials */
    var p;
    var c;
    if (parents.length==1) { // One arguments: glider on line, circle or curve
        p = parents[0];
        c = p.slideObject;
    } else if (parents.length==2) { // Two arguments: (point,line), (point,circle), (line,point) or (circle,point)
        if (JXG.isPoint(parents[0])) {
            p = parents[0];
            c = parents[1];
        } else if (JXG.isPoint(parents[1])) {
            c = parents[0];
            p = parents[1];
        } else {
            throw new Error("JSXGraph: Can't create normal with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
        }
    } else {
        throw new Error("JSXGraph: Can't create normal with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
    }

    if(c.elementClass==JXG.OBJECT_CLASS_LINE) {
        // return board.addNormal(c,p, attributes['id'], attributes['name']); // GEONExT-Style: problems with ideal point
        // If not needed, then board.addNormal and maybe board.algebra.perpendicular can be removed.

        // Homogeneous version:
        // orthogonal(l,p) = (F^\delta\cdot l)\times p
        return board.create('line', [
                    function(){ return c.stdform[1]*p.Y()-c.stdform[2]*p.X();},
                    function(){ return c.stdform[2]*p.Z();},
                    function(){ return -c.stdform[1]*p.Z();}
                    ], attributes );
    }
    else if(c.elementClass == JXG.OBJECT_CLASS_CIRCLE) {
        /*
        var Dg = function(t){ return -c.Radius()*Math.sin(t); };
        var Df = function(t){ return c.Radius()*Math.cos(t); };
        return board.create('line', [
                    function(){ return -p.X()*Dg(p.position)-p.Y()*Df(p.position);},
                    function(){ return Dg(p.position);},
                    function(){ return Df(p.position);}
                    ], attributes );
        */
        return board.create('line', [c.midpoint,p], attributes);
    } else if (c.elementClass == JXG.OBJECT_CLASS_CURVE) {
        if (c.curveType!='plot') {
            var g = c.X;
            var f = c.Y;
            return board.create('line', [
                    function(){ return -p.X()*board.D(g)(p.position)-p.Y()*board.D(f)(p.position);},
                    function(){ return board.D(g)(p.position);},
                    function(){ return board.D(f)(p.position);}
                    ], attributes );
        } else {                         // curveType 'plot'
            return board.create('line', [
                    function(){ var i=Math.floor(p.position);
                                var lbda = p.position-i;
                                if (i==c.numberPoints-1) {i--; lbda=1; }
                                if (i<0) return 1.0;
                                return (c.Y(i)+lbda*(c.Y(i+1)-c.Y(i)))*(c.Y(i)-c.Y(i+1))-(c.X(i)+lbda*(c.X(i+1)-c.X(i)))*(c.X(i+1)-c.X(i));},
                    function(){ var i=Math.floor(p.position);
                                if (i==c.numberPoints-1) i--;
                                if (i<0) return 0.0;
                                return c.X(i+1)-c.X(i);},
                    function(){ var i=Math.floor(p.position);
                                if (i==c.numberPoints-1) i--;
                                if (i<0) return 0.0;
                                return c.Y(i+1)-c.Y(i);}
                    ], attributes );
        }
    } else if (c.type == JXG.OBJECT_TYPE_TURTLE) {
            return board.create('line', [
                    function(){ var i=Math.floor(p.position);
                                var lbda = p.position-i;
                                var el,j;
                                for(j=0;j<c.objects.length;j++) {  // run through all curves of this turtle
                                    el = c.objects[j];
                                    if (el.type==JXG.OBJECT_TYPE_CURVE) {
                                        if (i<el.numberPoints) break;
                                        i-=el.numberPoints;
                                    }
                                }
                                if (i==el.numberPoints-1) { i--; lbda=1.0; }
                                if (i<0) return 1.0;
                                return (el.Y(i)+lbda*(el.Y(i+1)-el.Y(i)))*(el.Y(i)-el.Y(i+1))-(el.X(i)+lbda*(el.X(i+1)-el.X(i)))*(el.X(i+1)-el.X(i));},
                    function(){ var i=Math.floor(p.position);
                                var el,j;
                                for(j=0;j<c.objects.length;j++) {  // run through all curves of this turtle
                                    el = c.objects[j];
                                    if (el.type==JXG.OBJECT_TYPE_CURVE) {
                                        if (i<el.numberPoints) break;
                                        i-=el.numberPoints;
                                    }
                                }
                                if (i==el.numberPoints-1) i--;
                                if (i<0) return 0.0;
                                return el.X(i+1)-el.X(i);},
                    function(){ var i=Math.floor(p.position);
                                var el,j;
                                for(j=0;j<c.objects.length;j++) {  // run through all curves of this turtle
                                    el = c.objects[j];
                                    if (el.type==JXG.OBJECT_TYPE_CURVE) {
                                        if (i<el.numberPoints) break;
                                        i-=el.numberPoints;
                                    }
                                }
                                if (i==el.numberPoints-1) i--;
                                if (i<0) return 0.0;
                                return el.Y(i+1)-el.Y(i);}
                    ], attributes );
    }
    else {
        throw new Error("JSXGraph: Can't create normal with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
    }
};

/**
 * @class Provides a constructor for an angle bisector.
 * @pseudo
 * @description A bisector is a line which divides an angle into two equal angles. It is given by three points A, B, and C and divides the angle ABC into two
 * equal sized parts.
 * @constructor
 * @name Bisector
 * @type JXG.Line
 * @augments JXG.Line
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point_JXG.Point_JXG.Point} p1,p2,p3 The angle described by p3 will be divided into two equal angles.
 * @example
 * // Create a normal to a circle.
 * var p1 = board.create('point', [6.0, 4.0]);
 * var p2 = board.create('point', [3.0, 2.0]);
 * var p3 = board.create('point', [1.0, 7.0]);
 *
 * var bi1 = board.create('bisector', [p1, p2, p3]);
 * </pre><div id="0d58cea8-b06a-407c-b27c-0908f508f5a4" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var biex1_board = JXG.JSXGraph.initBoard('0d58cea8-b06a-407c-b27c-0908f508f5a4', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var biex1_p1 = biex1_board.create('point', [6.0, 4.0]);
 *   var biex1_p2 = biex1_board.create('point', [3.0, 2.0]);
 *   var biex1_p3 = biex1_board.create('point', [1.0, 7.0]);
 *   var biex1_bi1 = biex1_board.create('bisector', [biex1_p1, biex1_p2, biex1_p3]);
 * </script><pre>
 */
JXG.createBisector = function(board, parentArr, atts) {
    var p, l, cAtts, i;
    /* TODO bisector polynomials */
    if(parentArr[0].elementClass == JXG.OBJECT_CLASS_POINT && parentArr[1].elementClass == JXG.OBJECT_CLASS_POINT && parentArr[2].elementClass == JXG.OBJECT_CLASS_POINT) {

        cAtts = {name: '', id: null, fixed: true, visible: false};
        if(atts) {
            cAtts = JXG.cloneAndCopy(atts, cAtts);
        }

        // hidden and fixed helper
        p = board.create('point', [function () { return board.algebra.angleBisector(parentArr[0], parentArr[1], parentArr[2]); }], cAtts);

        for(i=0; i<3; i++)
            parentArr[i].addChild(p); // required for algorithm requiring dependencies between elements

        if(typeof atts['straightFirst'] == 'undefined')
            atts['straightFirst'] = false;
        if(typeof atts['straightLast'] == 'undefined')
            atts['straightLast'] = true;
        // no need to fire up the create stack because only attributes need to be set and they
        // will be set for l after returning.
        l = JXG.createLine(board, [parentArr[1], p], atts);
        return l;
    }
    else {
        throw new Error("JSXGraph: Can't create angle bisector with parent types '" + (typeof parentArr[0]) + "' and '" + (typeof parentArr[1]) + "'.");
    }
};

/**
 * TODO Is it possible to merge this with createBisector? --michael
 * The angular bisectors of two line [c1,a1,b1] and [c2,a2,b2] are determined by the equation:
 * (a1*x+b1*y+c1*z)/sqrt(a1^2+b1^2) = +/- (a2*x+b2*y+c2*z)/sqrt(a2^2+b2^2)
 * @private
 */
JXG.createAngularBisectorsOfTwoLines = function(board, parents, attributes) {
    var l1 = JXG.getReference(board,parents[0]),
        l2 = JXG.getReference(board,parents[1]),
        id1 = '',
        id2 = '',
        n1 = '',
        n2 = '',
        ret;

    attributes = JXG.checkAttributes(attributes,{});
    if (attributes['id']!=null) {
        if (JXG.isArray(attributes['id'])) {
            id1 = attributes['id'][0];
            id2 = attributes['id'][1];
        } else {
            id1 = attributes['id'];
            id2 = attributes['id'];
        }
    }
    if (attributes['name']!=null) {
        if (JXG.isArray(attributes['name'])) {
            n1 = attributes['name'][0];
            n2 = attributes['name'][1];
        } else {
            n1 = attributes['name'];
            n2 = attributes['name'];
        }
    }

    attributes['id'] = id1;
    attributes['name'] = n1;
    var g1 = board.create('line',[
        function(){
            var d1 = Math.sqrt(l1.stdform[1]*l1.stdform[1]+l1.stdform[2]*l1.stdform[2]);
            var d2 = Math.sqrt(l2.stdform[1]*l2.stdform[1]+l2.stdform[2]*l2.stdform[2]);
            return l1.stdform[0]/d1-l2.stdform[0]/d2;
        },
        function(){
            var d1 = Math.sqrt(l1.stdform[1]*l1.stdform[1]+l1.stdform[2]*l1.stdform[2]);
            var d2 = Math.sqrt(l2.stdform[1]*l2.stdform[1]+l2.stdform[2]*l2.stdform[2]);
            return l1.stdform[1]/d1-l2.stdform[1]/d2;
        },
        function(){
            var d1 = Math.sqrt(l1.stdform[1]*l1.stdform[1]+l1.stdform[2]*l1.stdform[2]);
            var d2 = Math.sqrt(l2.stdform[1]*l2.stdform[1]+l2.stdform[2]*l2.stdform[2]);
            return l1.stdform[2]/d1-l2.stdform[2]/d2;
        }
    ], attributes);
    attributes['id'] = id2;
    attributes['name'] = n2;
    var g2 = board.create('line',[
        function(){
            var d1 = Math.sqrt(l1.stdform[1]*l1.stdform[1]+l1.stdform[2]*l1.stdform[2]);
            var d2 = Math.sqrt(l2.stdform[1]*l2.stdform[1]+l2.stdform[2]*l2.stdform[2]);
            return l1.stdform[0]/d1+l2.stdform[0]/d2;
        },
        function(){
            var d1 = Math.sqrt(l1.stdform[1]*l1.stdform[1]+l1.stdform[2]*l1.stdform[2]);
            var d2 = Math.sqrt(l2.stdform[1]*l2.stdform[1]+l2.stdform[2]*l2.stdform[2]);
            return l1.stdform[1]/d1+l2.stdform[1]/d2;
        },
        function(){
            var d1 = Math.sqrt(l1.stdform[1]*l1.stdform[1]+l1.stdform[2]*l1.stdform[2]);
            var d2 = Math.sqrt(l2.stdform[1]*l2.stdform[1]+l2.stdform[2]*l2.stdform[2]);
            return l1.stdform[2]/d1+l2.stdform[2]/d2;
        }
    ], attributes);

    ret = [g1, g2];
    ret.lines = [g1, g2];
    ret.line1 = g1;
    ret.line2 = g2;

    ret.multipleElements = true;

    return ret;
};

/**
 * @class Constructs the midpoint of a {@link Circumcircle}. Like the circumcircle the circumcirclemidpoint
 * is constructed by providing three points.
 * @pseudo
 * @description A circumcircle midpoint is given by three points which are all lying on the circle with the
 * constructed circumcircle midpoint as the midpoint.
 * @constructor
 * @name Circumcirclemidpoint
 * @type JXG.Point
 * @augments JXG.Point
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point_JXG.Point_JXG.Point} p1,p2,p3 The constructed point is the midpoint of the circle determined
 * by p1, p2, and p3.
 * @example
 * var p1 = board.create('point', [0.0, 2.0]);
 * var p2 = board.create('point', [2.0, 1.0]);
 * var p3 = board.create('point', [3.0, 3.0]);
 *
 * var cc1 = board.create('circumcirclemidpoint', [p1, p2, p3]);
 * </pre><div id="e8a40f95-bf30-4eb4-88a8-f4d5495261fd" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var ccmex1_board = JXG.JSXGraph.initBoard('e8a40f95-bf30-4eb4-88a8-f4d5495261fd', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var ccmex1_p1 = ccmex1_board.create('point', [0.0, 2.0]);
 *   var ccmex1_p2 = ccmex1_board.create('point', [6.0, 1.0]);
 *   var ccmex1_p3 = ccmex1_board.create('point', [3.0, 7.0]);
 *   var ccmex1_cc1 = ccmex1_board.create('circumcirclemidpoint', [ccmex1_p1, ccmex1_p2, ccmex1_p3]);
 * </script><pre>
 */
JXG.createCircumcircleMidpoint = function(board, parentArr, atts) {
    var p, i;

    /* TODO circumcircle polynomials */

    if(parentArr[0].elementClass == JXG.OBJECT_CLASS_POINT && parentArr[1].elementClass == JXG.OBJECT_CLASS_POINT && parentArr[2].elementClass == JXG.OBJECT_CLASS_POINT) {
        atts['fixed'] = atts['fixed'] || true;
        p = JXG.createPoint(board, [function () { return board.algebra.circumcenterMidpoint(parentArr[0], parentArr[1], parentArr[2]); }], atts);

        for(i=0; i<3; i++)
            parentArr[i].addChild(p);

        return p;
    }
    else {
        throw new Error("JSXGraph: Can't create circumcircle midpoint with parent types '" + (typeof parentArr[0]) + "', '" + (typeof parentArr[1]) + "' and '" + (typeof parentArr[2]) + "'.");
    }
};

/**
 * @class Constructs two elements: a point and a circle. The circle is given by three points which lie on the circle,
 * the point is the midpoint of the circle.
 * @pseudo
 * @description A circumcircle is given by three points which are all lying on the circle.
 * @constructor
 * @name Circumcircle
 * @type array
 * @returns An array containing the midpoint in the first component and the circumcircle in the second component.
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point_JXG.Point_JXG.Point} p1,p2,p3 The constructed point is the midpoint of the circle determined
 * by p1, p2, and p3.
 * @example
 * var p1 = board.create('point', [0.0, 2.0]);
 * var p2 = board.create('point', [2.0, 1.0]);
 * var p3 = board.create('point', [3.0, 3.0]);
 *
 * var cc1 = board.create('circumcircle', [p1, p2, p3]);
 * </pre><div id="e65c9861-0bf0-402d-af57-3ab11962f5ac" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var ccex1_board = JXG.JSXGraph.initBoard('e65c9861-0bf0-402d-af57-3ab11962f5ac', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var ccex1_p1 = ccex1_board.create('point', [0.0, 2.0]);
 *   var ccex1_p2 = ccex1_board.create('point', [6.0, 1.0]);
 *   var ccex1_p3 = ccex1_board.create('point', [3.0, 7.0]);
 *   var ccex1_cc1 = ccex1_board.create('circumcircle', [ccex1_p1, ccex1_p2, ccex1_p3]);
 * </script><pre>
 */
JXG.createCircumcircle = function(board, parentArr, atts) {
    var p, c, cAtts, ret;

    cAtts = JXG.clone(atts);
    if(atts['name'] && JXG.isArray(atts['name'])) {
        cAtts['name'] = atts['name'][0];
        atts['name'] = atts['name'][1];
    }
    if(atts['id'] && JXG.isArray(atts['id'])) {
        cAtts['id'] = atts['id'][0];
        atts['id'] = atts['id'][1];
    }

    try {
        p = JXG.createCircumcircleMidpoint(board, parentArr, cAtts);
        c = JXG.createCircle(board, [p, parentArr[0]], atts);
    } catch(e) {
        throw new Error("JSXGraph: Can't create circumcircle with parent types '" + (typeof parentArr[0]) + "', '" + (typeof parentArr[1]) + "' and '" + (typeof parentArr[2]) + "'.");
    }

    ret = [p, c];

    ret.point = p;
    ret.circle = c;

    ret.multipleElements = true;

    return ret;
};

/**
 * @class This element is used to construct a reflected point.
 * @pseudo
 * @description A reflected point is given by a point and a line. It is determined by the reflection of the given point
 * against the given line.
 * @constructor
 * @name Reflection
 * @type JXG.Point
 * @augments JXG.Point
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point_JXG.Line} p,l The reflection point is the reflection of p against l.
 * @example
 * var p1 = board.create('point', [0.0, 4.0]);
 * var p2 = board.create('point', [6.0, 1.0]);
 * var l1 = board.create('line', [p1, p2]);
 * var p3 = board.create('point', [3.0, 3.0]);
 *
 * var rp1 = board.create('reflection', [p3, l1]);
 * </pre><div id="087a798e-a36a-4f52-a2b4-29a23a69393b" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var rpex1_board = JXG.JSXGraph.initBoard('087a798e-a36a-4f52-a2b4-29a23a69393b', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var rpex1_p1 = rpex1_board.create('point', [0.0, 4.0]);
 *   var rpex1_p2 = rpex1_board.create('point', [6.0, 1.0]);
 *   var rpex1_l1 = rpex1_board.create('line', [rpex1_p1, rpex1_p2]);
 *   var rpex1_p3 = rpex1_board.create('point', [3.0, 3.0]);
 *   var rpex1_rp1 = rpex1_board.create('reflection', [rpex1_p3, rpex1_l1]);
 * </script><pre>
 */
JXG.createReflection = function(board, parentArr, atts) {
    var l, p, r;

    /* TODO reflection polynomials */
    if(parentArr[0].elementClass == JXG.OBJECT_CLASS_POINT && parentArr[1].elementClass == JXG.OBJECT_CLASS_LINE) {
        p = parentArr[0];
        l = parentArr[1];
    }
    else if(parentArr[1].elementClass == JXG.OBJECT_CLASS_POINT && parentArr[0].elementClass == JXG.OBJECT_CLASS_LINE) {
        p = parentArr[1];
        l = parentArr[0];
    }
    else {
        throw new Error("JSXGraph: Can't create reflection point with parent types '" + (typeof parentArr[0]) + "' and '" + (typeof parentArr[1]) + "'.");
    }

    // force a fixed point
    atts['fixed'] = true;
    r = JXG.createPoint(board, [function () { return board.algebra.reflection(l, p); }], atts);
    p.addChild(r);
    l.addChild(r);

    r.update();

    return r;
};

// here we have to continue with replacing board.add* stuff

/**
 * @class A mirror point will be constructed.
 * @pseudo
 * @description A mirror point is determined by the reflection of a given point against another given point.
 * @constructor
 * @name Mirrorpoint
 * @type JXG.Point
 * @augments JXG.Point
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point_JXG.Point} p1,p2 The constructed point is the reflection of p2 against p1.
 * @example
 * var p1 = board.create('point', [3.0, 3.0]);
 * var p2 = board.create('point', [6.0, 1.0]);
 *
 * var mp1 = board.create('mirrorpoint', [p1, p2]);
 * </pre><div id="7eb2a814-6c4b-4caa-8cfa-4183a948d25b" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var mpex1_board = JXG.JSXGraph.initBoard('7eb2a814-6c4b-4caa-8cfa-4183a948d25b', {boundingbox: [-1, 9, 9, -1], axis: true, showcopyright: false, shownavigation: false});
 *   var mpex1_p1 = mpex1_board.create('point', [3.0, 3.0]);
 *   var mpex1_p2 = mpex1_board.create('point', [6.0, 1.0]);
 *   var mpex1_mp1 = mpex1_board.create('mirrorpoint', [mpex1_p1, mpex1_p2]);
 * </script><pre>
 */
JXG.createMirrorPoint = function(board, parentArr, atts) {
    var p;

    /* TODO mirror polynomials */
    if(JXG.isPoint(parentArr[0]) && JXG.isPoint(parentArr[1])) {
        atts['fixed'] = atts['fixed'] || true;
        p = JXG.createPoint(board, [function () { return board.algebra.rotation(parentArr[0], parentArr[1], Math.PI); }], atts);

        for(i=0; i<2; i++)
            parentArr[i].addChild(p);
    }
    else {
        throw new Error("JSXGraph: Can't create mirror point with parent types '" + (typeof parentArr[0]) + "' and '" + (typeof parentArr[1]) + "'.");
    }

    p.update();

    return p;
};

/**
 * @class This element is used to visualize the integral of a given curve over a given interval.
 * @pseudo
 * @description The Integral element is used to visualize the area under a given curve over a given interval
 * and to calculate the area's value. For that a polygon and gliders are used. The polygon displays the area,
 * the gliders are used to change the interval dynamically.
 * @constructor
 * @name Integral
 * @type JXG.Polygon
 * @augments JXG.Polygon
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {array_JXG.Curve} p,l The constructed point is the orthogonal projection of p onto l.
 * @example
 * var c1 = board.create('functiongraph', [function (t) { return t*t*t; }]);
 * var i1 = board.create('integral', [[-1.0, 4.0], c1]);
 * </pre><div id="d45d7188-6624-4d6e-bebb-1efa2a305c8a" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var intex1_board = JXG.JSXGraph.initBoard('d45d7188-6624-4d6e-bebb-1efa2a305c8a', {boundingbox: [-5, 5, 5, -5], axis: true, showcopyright: false, shownavigation: false});
 *   var intex1_c1 = intex1_board.create('functiongraph', [function (t) { return t*t*t; }]);
 *   var intex1_i1 = intex1_board.create('integral', [[-2.0, 2.0], intex1_c1]);
 * </script><pre>
 */
JXG.createIntegral = function(board, parents, attributes) {
    var interval, curve, attribs = {},
        start = 0, end = 0,
        pa_on_curve, pa_on_axis, pb_on_curve, pb_on_axis,
        Int, t, p;

    if(!JXG.isArray(attributes['id']) || (attributes['id'].length != 5)) {
        attributes['id'] = ['','','','',''];
    }
    if(!JXG.isArray(attributes['name']) || (attributes['name'].length != 5)) {
       attributes['name'] = ['','','','',''];
    }

    if(JXG.isArray(parents[0]) && parents[1].type == JXG.OBJECT_TYPE_CURVE) {
        interval = parents[0];
        curve = parents[1];
    } else if(JXG.isArray(parents[1]) && parents[0].type == JXG.OBJECT_TYPE_CURVE) {
        interval = parents[1];
        curve = parents[0];
    } else {
        throw new Error("JSXGraph: Can't create integral with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "'.");
    }

    if( (typeof attributes != 'undefined') && (attributes != null))
        attribs = JXG.cloneAndCopy(attributes, {name: attributes.name[0], id: attributes.id[0]});

    // Correct the interval if necessary
    if(interval[0] > curve.points[0].usrCoords[1])
        start = interval[0];
    else
        start = curve.points[0].usrCoords[1];

    if(interval[1] < curve.points[curve.points.length-1].usrCoords[1])
        end = interval[1];
    else
        end = curve.points[curve.points.length-1].usrCoords[1];

    pa_on_curve = board.create('glider', [start, curve.yterm(start), curve], attribs);

    attribs.name = attributes.name[1];
    attribs.id = attributes.id[1];
    attribs.visible = false;
    pa_on_axis = board.create('point', [function () { return pa_on_curve.X(); }, 0], attribs);

    pa_on_curve.addChild(pa_on_axis);

    attribs.name = attributes.name[2];
    attribs.id = attributes.id[2];
    attribs.visible = attributes.visible || true;
    pb_on_curve = board.create('glider', [end, curve.yterm(end), curve], attribs);

    attribs.name = attributes.name[3];
    attribs.id = attributes.id[3];
    attribs.visible = false;
    pb_on_axis = board.create('point', [function () { return pb_on_curve.X(); }, 0], attribs);

    pb_on_curve.addChild(pb_on_axis);

    Int = JXG.Math.Numerics.I([start, end], curve.yterm);
    t = board.create('text', [
        function () { return pb_on_curve.X() + 0.2; },
        function () { return pb_on_curve.Y() - 0.8; },
        function () {
                var Int = JXG.Math.Numerics.I([pa_on_axis.X(), pb_on_axis.X()], curve.yterm);
                return '&int; = ' + (Int).toFixed(4);
            }
        ],{labelColor: attributes['labelColor']});

    attribs.name = attributes.name[4];
    attribs.id = attributes.id[4];
    attribs.visible = attributes.visible || true;
    attribs.fillColor = attribs.fillColor || board.options.polygon.fillColor;
    attribs.highlightFillColor = attribs.highlightFillColor || board.options.polygon.highlightFillColor;
    attribs.fillOpacity = attribs.fillOpacity || board.options.polygon.fillOpacity;
    attribs.highlightFillOpacity = attribs.highlightFillOpacity || board.options.polygon.highlightFillOpacity;
    attribs.strokeWidth = 0;
    attribs.strokeOpacity = 0;

    p = board.create('curve', [[0],[0]], attribs);
    p.updateDataArray = function() {
        var x = [pa_on_axis.coords.usrCoords[1], pa_on_curve.coords.usrCoords[1]],
            y = [pa_on_axis.coords.usrCoords[2], pa_on_curve.coords.usrCoords[2]],
            i;

        for(i=0; i < curve.numberPoints; i++) {
            if( (pa_on_axis.X() <= curve.points[i].usrCoords[1]) && (curve.points[i].usrCoords[1] <= pb_on_axis.X()) ) {
                x.push(curve.points[i].usrCoords[1]);
                y.push(curve.points[i].usrCoords[2]);
            }
        }
        x.push(pb_on_curve.coords.usrCoords[1]);
        y.push(pb_on_curve.coords.usrCoords[2]);
        x.push(pb_on_axis.coords.usrCoords[1]);
        y.push(pb_on_axis.coords.usrCoords[2]);

        x.push(pa_on_axis.coords.usrCoords[1]); // close the curve
        y.push(pa_on_axis.coords.usrCoords[2]);

        this.dataX = x;
        this.dataY = y;
    }
    pa_on_curve.addChild(p);
    pb_on_curve.addChild(p);
    pa_on_curve.addChild(t);
    pb_on_curve.addChild(t);

    return p;//[pa_on_axis, pb_on_axis, p, t];

};

/**
 * @class This element is used to visualize the locus of a given dependent point.
 * @pseudo
 * @description The locus element is used to visualize the curve a given point describes.
 * @constructor
 * @name Locus
 * @type JXG.Curve
 * @augments JXG.Curve
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {JXG.Point} p The constructed curve is the geometric locus of the given point.
 * @example
 *  // This examples needs JXG.Server up and running, otherwise it won't work.
 *  p1 = board.create('point', [0, 0]);
 *  p2 = board.create('point', [6, -1]);
 *  c1 = board.create('circle', [p1, 2]);
 *  c2 = board.create('circle', [p2, 1.5]);
 *  g1 = board.create('glider', [6, 3, c1]);
 *  c3 = board.create('circle', [g1, 4]);
 *  g2 = board.create('intersection', [c2,c3,0]);
 *  m1 = board.create('midpoint', [g1,g2]);
 *  loc = board.create('locus', [m1], {strokeColor: 'red'});
 * </pre><div id="d45d7188-6624-4d6e-bebb-1efa2a305c8a" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *  lcex_board = JXG.JSXGraph.initBoard('jxgbox', {boundingbox:[-4, 6, 10, -6], axis: true, grid: false, keepaspectratio: true});
 *  lcex_p1 = lcex_board.create('point', [0, 0]);
 *  lcex_p2 = lcex_board.create('point', [6, -1]);
 *  lcex_c1 = lcex_board.create('circle', [lcex_p1, 2]);
 *  lcex_c2 = lcex_board.create('circle', [lcex_p2, 1.5]);
 *  lcex_g1 = lcex_board.create('glider', [6, 3, lcex_c1]);
 *  lcex_c3 = lcex_board.create('circle', [lcex_g1, 4]);
 *  lcex_g2 = lcex_board.create('intersection', [lcex_c2,lcex_c3,0]);
 *  lcex_m1 = lcex_board.create('midpoint', [lcex_g1,lcex_g2]);
 *  lcex_loc = board.create('locus', [lcex_m1], {strokeColor: 'red'});
 * </script><pre>
 */
JXG.createLocus = function(board, parents, attributes) {
    var c, p;

    if(JXG.isArray(parents) && parents.length == 1 && parents[0].elementClass == JXG.OBJECT_CLASS_POINT) {
        p = parents[0];
    } else {
        throw new Error("JSXGraph: Can't create locus with parent of type other than point.");
    }

    c = board.create('curve', [[null], [null]], attributes);
    c.dontCallServer = false;

    c.updateDataArray = function () {
        cb = function(x, y, eq) {
            c.dataX = x;
            c.dataY = y;
            board.update();
        };

        if(board.mode == board.BOARD_MODE_NONE && !this.dontCallServer) {
            JXG.Math.Symbolic.geometricLocusByGroebnerBase(board, p, cb);
            // don't bother the server on the next update, because it's fired
            // to plot the datapoints received by the server.
            this.dontCallServer = true;
        } else {
            this.dontCallServer = false;
        }
    };
    return c;
};

JXG.JSXGraph.registerElement('arrowparallel', JXG.createArrowParallel);
JXG.JSXGraph.registerElement('bisector', JXG.createBisector);
JXG.JSXGraph.registerElement('bisectorlines', JXG.createAngularBisectorsOfTwoLines);
JXG.JSXGraph.registerElement('circumcircle', JXG.createCircumcircle);
JXG.JSXGraph.registerElement('circumcirclemidpoint', JXG.createCircumcircleMidpoint);
JXG.JSXGraph.registerElement('integral', JXG.createIntegral);
JXG.JSXGraph.registerElement('midpoint', JXG.createMidpoint);
JXG.JSXGraph.registerElement('mirrorpoint', JXG.createMirrorPoint);
JXG.JSXGraph.registerElement('normal', JXG.createNormal);
JXG.JSXGraph.registerElement('parallel', JXG.createParallel);
JXG.JSXGraph.registerElement('parallelpoint', JXG.createParallelPoint);
JXG.JSXGraph.registerElement('perpendicular', JXG.createPerpendicular);
JXG.JSXGraph.registerElement('perpendicularpoint', JXG.createPerpendicularPoint);
JXG.JSXGraph.registerElement('reflection', JXG.createReflection);
JXG.JSXGraph.registerElement('locus', JXG.createLocus);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview In this file the Text element is defined.
 */


/**
 * Construct and handle texts.
 * @class Text: On creation the GEONExT syntax
 * of <value>-terms 
 * are converted into JavaScript syntax.
 * The coordinates can be relative to the coordinates of an element "element".
 * @constructor
 * @return A new geometry element Text
 */
JXG.Text = function (board, contentStr, element, coords, id, name, digits, isLabel, display, layer) {
    this.constructor();

    this.type = JXG.OBJECT_TYPE_TEXT;
    this.elementClass = JXG.OBJECT_CLASS_OTHER;                

    this.init(board, id, name);

    this.contentStr = contentStr;
    this.plaintextStr = '';

    /**
     * Set the display layer.
     */
    if (layer == null) layer = board.options.layer['text'];
    this.layer = layer;

    /**
     * There is choice between 'html' and 'internal'
     * 'internal' is the text element of SVG and the textpath element 
     * of VML.
     */
    this.display = display || 'html'; 
    
    if((typeof isLabel != 'undefined') && (isLabel != null)) {
        this.isLabel = isLabel;
    }
    else {
        this.isLabel = false;
    }
    
    /**
     * The text color of the given text.
     * @type {string}
     * @name JXG.Text#strokeColor
     */
    this.visProp['strokeColor'] = this.board.options.text.strokeColor;
    /**
     * The text opacity of the given text.
     * @type {string}
     * @name JXG.Text#strokeOpacity
     */
     /**
     * The font size of the given text.
     * @type {string}
     * @name JXG.Text#fontSize
     * @default {@link JXG.Options.fontSize}
     */

    this.visProp['visible'] = true;
    //this.show = true; // noch noetig? BV

    if (digits!=null) {
        this.digits = digits;
    } else {
        this.digits = 2;
    }

    /**
     * Coordinates of the text.
     * @ private
     * @type JXG.Coords
     */
    if ((this.element = this.board.objects[element])){
        var anchor;
        this.relativeCoords = new JXG.Coords(JXG.COORDS_BY_USER, [parseFloat(coords[0]),parseFloat(coords[1])],this.board);     
        if(!this.isLabel) {
            anchor = this.element.getTextAnchor();
        }
        else {
            anchor = this.element.getLabelAnchor();
        }      
        this.element.addChild(this);
        this.coords = new JXG.Coords(JXG.COORDS_BY_USER, [this.relativeCoords.usrCoords[1]+anchor.usrCoords[1],this.relativeCoords.usrCoords[2]+anchor.usrCoords[2]], this.board);
    } else {
        this.X = JXG.createFunction(coords[0],this.board,'');
        this.Y = JXG.createFunction(coords[1],this.board,'');
        this.coords = new JXG.Coords(JXG.COORDS_BY_USER, [this.X(),this.Y()], this.board);
        var fs = 'this.coords.setCoordinates(JXG.COORDS_BY_USER,[this.X(),this.Y()]);';
        this.updateCoords = new Function('',fs);
    }

    if (typeof this.contentStr=='function') {
        this.updateText = function() { this.plaintextStr = this.contentStr(); };
    } else {
        var plaintext;
        if (typeof this.contentStr=='number') {
            plaintext = (this.contentStr).toFixed(this.digits);  
        } else {
            if (this.board.options.text.useASCIIMathML) {
                plaintext = "'"+this.contentStr+"'";              // Convert via ASCIIMathML
            } else {
                plaintext = this.generateTerm(this.contentStr);   // Converts GEONExT syntax into JavaScript string
            }
        }
        this.updateText = new Function('this.plaintextStr = ' + plaintext + ';');
    }
    //this.updateText();                    // First evaluation of the string    
    if(!this.isLabel) {
        this.id = this.board.addText(this);
    }
    if (typeof this.contentStr=='string') {
        this.notifyParents(this.contentStr);
    }
};
JXG.Text.prototype = new JXG.GeometryElement();

/**
 * @private
 * Empty function (for the moment). It is needed for highlighting
 * @param {int} x
 * @param {int} y Find closest point on the text to (xy)
 * @return Always returns false
 */
JXG.Text.prototype.hasPoint = function (x,y) {
    return false;
};

/**
 * Overwrite the text.
 * @param {string,function} str
 * @return {object} reference to the text object.
 */
JXG.Text.prototype.setText = function(text) {
    var plaintext;
    if (typeof text=='number') {
        plaintext = (text).toFixed(this.digits);  
    } else {
        plaintext = this.generateTerm(text);   // Converts GEONExT syntax into JavaScript string
    }
    this.updateText = new Function('this.plaintextStr = ' + plaintext + ';');
    this.updateText();
    return this;
};

/**
 * Set the text to new, fixed coordinates.
 * @param {number} x
 * @param {number} y
 * @return {object} reference to the text object.
 */
JXG.Text.prototype.setCoords = function (x,y) {
    this.X = function() { return x; };
    this.Y = function() { return y; };
    this.coords = new JXG.Coords(JXG.COORDS_BY_USER, [x,y], this.board);
    return this;
};

/**
 * Evaluates the text.
 * Then, the update function of the renderer
 * is called. 
 */
JXG.Text.prototype.update = function () {
    if (this.needsUpdate) {
        if (this.relativeCoords){
            var anchor;
            if(!this.isLabel) {
                anchor = this.element.getTextAnchor();
            }
            else {
                anchor = this.element.getLabelAnchor();
            }
            this.coords.setCoordinates(JXG.COORDS_BY_USER, [this.relativeCoords.usrCoords[1]+anchor.usrCoords[1],this.relativeCoords.usrCoords[2]+anchor.usrCoords[2]]);
        } else {
            this.updateCoords();
        }
        this.updateText();
    }   
    return this;
};

/**
 * The update function of the renderer
 * is called. 
 * @private
 */
JXG.Text.prototype.updateRenderer = function () {
    if (this.needsUpdate) {
        this.board.renderer.updateText(this);
        this.needsUpdate = false;
    }
    return this;
};

/**
 * Converts the GEONExT syntax of the <value> terms into JavaScript.
 * Also, all Objects whose name appears in the term are searched and
 * the text is added as child to these objects.
 * @private
 * @see Algebra
 * @see #geonext2JS.
 */
JXG.Text.prototype.generateTerm = function (contentStr) {
    var res = null;
    var elements = this.board.elementsByName;
    var plaintext = '""';
    contentStr = contentStr.replace(/\r/g,''); 
    contentStr = contentStr.replace(/\n/g,''); 
    contentStr = contentStr.replace(/\"/g,'\\"'); 
    contentStr = contentStr.replace(/\'/g,"\\'"); 
    contentStr = contentStr.replace(/&amp;arc;/g,'&ang;'); 
    contentStr = contentStr.replace(/<arc\s*\/>/g,'&ang;'); 
    contentStr = contentStr.replace(/<sqrt\s*\/>/g,'&radic;'); 

    // Convert GEONExT syntax into  JavaScript syntax
    var i;
    //var i = contentStr.indexOf('<mp>');
    //contentStr = contentStr.slice(i+4);
    //i = contentStr.indexOf('</mp>');
    //contentStr = contentStr.slice(0,i);

    i = contentStr.indexOf('<value>');
    var j = contentStr.indexOf('</value>');
    if (i>=0) {
        while (i>=0) {
            plaintext += ' + "'+ this.board.algebra.replaceSub(this.board.algebra.replaceSup(contentStr.slice(0,i))) + '"';
            var term = contentStr.slice(i+7,j);
            var res = this.board.algebra.geonext2JS(term); 
            res = res.replace(/\\"/g,'"');
            res = res.replace(/\\'/g,"'");
            if (res.indexOf('toFixed')<0) {  // GEONExT-Hack: apply rounding once only.  
                plaintext += '+('+ res + ').toFixed('+(this.digits)+')';
            } else {
                plaintext += '+('+ res + ')';
            }
            contentStr = contentStr.slice(j+8);
            i = contentStr.indexOf('<value>');
            j = contentStr.indexOf('</value>');
        }
    } //else {
    plaintext += ' + "' + this.board.algebra.replaceSub(this.board.algebra.replaceSup(contentStr)) + '"';
    //}
    plaintext = plaintext.replace(/<overline>/g,'<span style=text-decoration:overline>');
    plaintext = plaintext.replace(/<\/overline>/g,'</span>');
    plaintext = plaintext.replace(/<arrow>/g,'<span style=text-decoration:overline>');
    plaintext = plaintext.replace(/<\/arrow>/g,'</span>');

/*    i = plaintext.indexOf('<name>');
    j = plaintext.indexOf('</name>');
    while (i>=0) {
        var head = plaintext.slice(0,i+6);
        var mid = plaintext.slice(i+6,j);
        var tail = plaintext.slice(j);
        mid = this.board.algebra.replaceSub(this.board.algebra.replaceSup(mid));
        plaintext = head + mid + tail;
        i = plaintext.indexOf('<name>',i+7);
        j = plaintext.indexOf('</name>',i+7);
    }
*/
    plaintext = plaintext.replace(/&amp;/g,'&'); // This should replace &amp;pi; by &pi;
    return plaintext;
};

/**
 * Finds dependencies in a given term and notifies the parents by adding the
 * dependent object to the found objects child elements.
 * @param {String} term String containing dependencies for the given object.
 * @private
 */
JXG.Text.prototype.notifyParents = function (contentStr) {
    var res = null;
    var elements = this.board.elementsByName;

    do {
        var search = /<value>([\w\s\*\/\^\-\+\(\)\[\],<>=!]+)<\/value>/;
        res = search.exec(contentStr);
        if (res!=null) {
            this.board.algebra.findDependencies(this,res[1]);
            contentStr = contentStr.substr(res.index);
            contentStr = contentStr.replace(search,'');
        }
    } while (res!=null);
    return this;
};

/**
 * @class This element is used to provide a constructor for text, which is just a wrapper for element {@link Text}. 
 * @pseudo
 * @description
 * @name Text
 * @augments JXG.GeometryElement
 * @constructor
 * @type JXG.Text
 *
 * @param {number,function_number,function_String,function} x,y,str Parent elements for text elements.
 *                     <p>
 *                     x and y are the coordinates of the lower left corner of the text box. The position of the text is fixed, 
 *                     x and y are numbers. The position is variable if x or y are functions.
 *                     <p>
 *                     The text to display may be given as string or as function returning a string.
 *
 * There is the attribute 'display' which takes the values 'html' or 'internal'. In case of 'html' a HTML division tag is created to display
 * the text. In this case it is also possible to use ASCIIMathML. Incase of 'internal', a SVG or VML text element is used to display the text.
 * @see JXG.Text
 * @example
 * // Create a fixed text at position [0,1].
 *   var t1 = board.create('text',[0,1,"Hello World"]); 
 * </pre><div id="896013aa-f24e-4e83-ad50-7bc7df23f6b7" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var t1_board = JXG.JSXGraph.initBoard('896013aa-f24e-4e83-ad50-7bc7df23f6b7', {boundingbox: [-3, 6, 5, -3], axis: true, showcopyright: false, shownavigation: false});
 *   var t1 = t1_board.create('text',[0,1,"Hello World"]);
 * </script><pre>
 * @example
 * // Create a variable text at a variable position.
 *   var s = board.create('slider',[[0,4],[3,4],[-2,0,2]]);
 *   var graph = board.create('text', 
 *                        [function(x){ return s.Value();}, 1,
 *                         function(){return "The value of s is"+s.Value().toFixed(2);}
 *                        ]
 *                     );
 * </pre><div id="5441da79-a48d-48e8-9e53-75594c384a1c" style="width: 300px; height: 300px;"></div>
 * <script type="text/javascript">
 *   var t2_board = JXG.JSXGraph.initBoard('5441da79-a48d-48e8-9e53-75594c384a1c', {boundingbox: [-3, 6, 5, -3], axis: true, showcopyright: false, shownavigation: false});
 *   var s = t2_board.create('slider',[[0,4],[3,4],[-2,0,2]]);
 *   var t2 = t2_board.create('text',[function(x){ return s.Value();}, 1, function(){return "The value of s is "+s.Value().toFixed(2);}]);
 * </script><pre>
 */
JXG.createText = function(board, parentArr, atts) {
    atts = JXG.checkAttributes(atts,{layer:null,display:board.options.text.defaultDisplay});  // 'html' or 'internal'
    return new JXG.Text(board, parentArr[parentArr.length-1], null, parentArr, atts['id'], atts['name'], atts['digits'], false, atts['display'],atts['layer']);
};

JXG.JSXGraph.registerElement('text', JXG.createText);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.

*/

/**
 * @fileoverview In this file the geometry element Image is defined.
 * @author graphjs
 * @version 0.1
 */

/**
 * Construct and handle images
 * @class Image:
 * It inherits from @see GeometryElement.
 * @constructor
 * @return A new geometry element Image
 */
JXG.Image = function (board, url, coordinates, size, layer, id, name, el) {
    //this.constructor();
    this.type = JXG.OBJECT_TYPE_IMAGE;
    this.elementClass = JXG.OBJECT_CLASS_OTHER;                
    this.transformations = [];

    this.init(board, id, name);
    this.coords = new JXG.Coords(JXG.COORDS_BY_USER, coordinates, this.board);
    this.initialCoords = new JXG.Coords(JXG.COORDS_BY_USER, coordinates, this.board);
    this.size = [size[0]*board.stretchX,size[1]*board.stretchY];
    //this.imageBase64String = url; //imageBase64String;
    this.url = url;
    /**
     * Set the display layer.
     */
    if (layer == null) layer = board.options.layer['image'];
    this.layer = layer;
    this.parent = el;
    this.visProp['visible'] = true;
    //this.show = true; // noch noetig? BV
    this.id = this.board.addImage(this);
};

JXG.Image.prototype = new JXG.GeometryElement;

/**
 * Empty function (for the moment). It is needed for highlighting, a feature not used for images right now.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @return Always returns false
 */
JXG.Image.prototype.hasPoint = function (x,y) {
    return false;
};

/**
 * Send an update request to the renderer.
 */
JXG.Image.prototype.updateRenderer = function () {
    this.updateTransform();
    this.board.renderer.updateImage(this);
};

JXG.Image.prototype.updateTransform = function () {
    if (this.transformations.length==0) {
        return;
    }
    for (var i=0;i<this.transformations.length;i++) {
        this.transformations[i].update();
    }
};

JXG.Image.prototype.addTransform = function (transform) {
    if (JXG.isArray(transform)) {
        for (var i=0;i<transform.length;i++) {
            this.transformations.push(transform[i]);
        }
    } else {
        this.transformations.push(transform);
    }
};

/**
 * @class Displays an image. 
 * @pseudo
 * @description Shows an imgae. The image can be supplied as an URL or an base64 encoded inline image
 * like "data:image/png;base64, /9j/4AAQSkZJRgA...".
 * @constructor
 * @name Image
 * @type JXG.Image
 * @throws {Exception} If the element cannot be constructed with the given parent objects an exception is thrown.
 * @param {String_Array_Array} url, [position of the top left vertice], [width,height] 
 * @example
 * var im = board.create('image', ['http://geonext.uni-bayreuth.de/fileadmin/geonext/design/images/logo.gif', [-3,1],[5,5]]);
 *
 * </pre><div id="9850cda0-7ea0-4750-981c-68bacf9cca57" style="width: 400px; height: 400px;"></div>
 * <script type="text/javascript">
 *   var image_board = JXG.JSXGraph.initBoard('9850cda0-7ea0-4750-981c-68bacf9cca57', {boundingbox: [-4, 4, 4, -4], axis: false, showcopyright: false, shownavigation: false});
 *   var image_im = image_board.create('image', ['http://geonext.uni-bayreuth.de/fileadmin/geonext/design/images/logo.gif', [-3,1],[5,5]]);
 * </script><pre>
 */
JXG.createImage = function(board, parents, atts) {
    var url;
    if (atts==null) {
        atts = {};
    } else if (atts['imageString']!=null) {
        url = atts['imageString'];
    }
    if (typeof atts['layer'] == 'undefined') {
        atts['layer'] = null;
    }
    return new JXG.Image(board, parents[0], parents[1], parents[2], atts['layer'], false, false, undefined);
};

JXG.JSXGraph.registerElement('image', JXG.createImage);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview The geometry object Line is defined in this file. Line stores all
 * style and functional properties that are required to draw and move a line on
 * a board.
 */

/**
 * The Slider class is....
 * Slider (Schieberegler)
 * input: 3 arrays:
 * [x0,y0],[x1,y1],[min,start,max]
 * The slider is line from [x0,y0] to [x1,y1].
 * The position [x0,y0]  corresponds to the value "min",
 * [x1,y1] corresponds to the value max.
 * Initally, the slider is at position [x0,y0] + ([x1,y1]-[x0,y0])*start/(max-min)
 * The return value is an invisible point, whos X() or Y() value
 * returns the position between max and min,
 * Further, there is a method Value() returning the same value.
 * @class Creates a new basic slider object. Do not use this constructor to create a slider. Use {@link JXG.Board#create} with
 * type {@link Line}, {@link Arrow}, or {@link Axis} instead.  
 * @constructor
 * @augments JXG.GeometryElement
 * @param {String,JXG.Board} board The board the new line is drawn on.
 * @param {Point} p1 Startpoint of the line.
 * @param {Point} p2 Endpoint of the line.
 * @param {String} id Unique identifier for this object. If null or an empty string is given,
 * an unique id will be generated by Board
 * @param {String} name Not necessarily unique name. If null or an
 * empty string is given, an unique name will be generated.
 * @see JXG.Board#generateName
 */
JXG.createSlider = function(board, parentArr, atts) {
    var pos0, pos1, smin, start, smax, sdiff, p1, p2, l1, ticks, ti, startx, starty, p3, l2, n, t,
        snapWidth;
        
    pos0 = parentArr[0];
    pos1 = parentArr[1];
    smin = parentArr[2][0];
    start = parentArr[2][1];
    smax = parentArr[2][2];
    sdiff = smax -smin;
    
    atts = JXG.checkAttributes(atts,{strokeColor:'#000000', fillColor:'#ffffff'});

    p1 = board.create('point', pos0, {visible:false, fixed:true, name:'',withLabel:false}); 
    p2 = board.create('point', pos1, {visible:false, fixed:true, name:'',withLabel:false}); 
    l1 = board.create('segment', [p1,p2], 
                {strokewidth:1, 
                name:'',
                withLabel:false,
                strokeColor:atts['strokeColor']});
    ticks  = 2;
    ti = board.create('ticks', [l1, p2.Dist(p1)/ticks],
                {insertTicks:true, minorTicks:0, drawLabels:false, drawZero:true}); 

    p1.needsRegularUpdate = false;
    p2.needsRegularUpdate = false;
    l1.needsRegularUpdate = false;
    
    startx = pos0[0]+(pos1[0]-pos0[0])*(start-smin)/(smax-smin);
    starty = pos0[1]+(pos1[1]-pos0[1])*(start-smin)/(smax-smin);

    if (atts['snapWidth']!=null) snapWidth = atts['snapWidth'];
    if (atts['snapwidth']!=null) snapWidth = atts['snapwidth'];
    
    p3 = board.create('glider', [startx,starty,l1],
                {style:6,strokeColor:atts['strokeColor'],
                 fillColor:atts['fillColor'],
                 showInfobox:false,name:atts['name'], withLabel:false,
                 snapWidth:snapWidth});
    
    l2 = board.create('line', [p1,p3], 
                {straightFirst:false, 
                 straightLast:false, strokewidth:3, 
                 strokeColor:atts['strokeColor'],
                 name:'',
                 withLabel:false}); 
                 
    //p3.Value = function() { return this.position*(smax - smin)+smin; };
    //p3.type = JXG.OBJECT_TYPE_SLIDER;
    p3.Value = function() { return this.position*sdiff+smin; };
    p3._smax = smax;
    p3._smin = smin;

    if (typeof atts['withLabel']=='undefined' || atts['withLabel']==true) {
        if (atts['name'] && atts['name']!='') {
            n = atts['name'] + ' = ';
        } else {
            n = '';
        }
        t = board.create('text', [((pos1[0]-pos0[0])*.05+pos1[0]), 
                                     ((pos1[1]-pos0[1])*.05+pos1[1]), 
                                     function(){return n+(p3.Value()).toFixed(2);}],
                                     {name:''}); 
    }                                     
    return p3;
};    

JXG.JSXGraph.registerElement('slider', JXG.createSlider);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * Chart plotting
 * input array: 
 * 
 **/
 
JXG.Chart = function(board, parents, attributes) {
    this.constructor();
    if (parents.length==0) { return; }  // No input data in parentArr
    
    /**
     * Contains lpointers to the various displays.
     */
    this.elements = [];
    
    var id = attributes['id'] || '';
    var name = attributes['name'] || '';
    this.init(board, id, name);
    
    var x,y,i;
    if (parents.length>0 && (typeof parents[0]=='number')) { // parents looks like [a,b,c,..]
                                                                // x has to be filled
        y = parents;
        x = [];
        for (i=0;i<y.length;i++) {
            x[i] = i+1;
        }
    } else {
        if (parents.length==1) { // parents looks like [[a,b,c,..]]
                                // x has to be filled
            y = parents[0];
            x = [];
            var len;
            if (JXG.isFunction(y)) {
                len = y().length;
            } else {
                len = y.length;
            }
            for (i=0;i<len;i++) {
                x[i] = i+1;
            }
        }
        if (parents.length==2) { // parents looks like [[x0,x1,x2,...],[y1,y2,y3,...]]
            x = parents[0];
            y = parents[1];
        }
    }
    if (attributes==undefined) attributes = {};
    var style = attributes['chartStyle'] || 'line';
    style = style.replace(/ /g,'');
    style = style.split(',');
    var c;
    for (i=0;i<style.length;i++) {
        switch (style[i]) {
            case 'bar':
                c = this.drawBar(board,[x,y],attributes);
                break;
            case 'line':
                c = this.drawLine(board, [x, y], attributes);
                break;
            case 'fit':
                c = this.drawFit(board, [x, y], attributes);
                break;
            case 'spline':
                c = this.drawSpline(board, [x, y], attributes);
                break;
            case 'pie':
                c = this.drawPie(board,[y],attributes);
                break;
            case 'point':
                c = this.drawPoints(board,[x,y],attributes);
                break;
        };
        this.elements.push(c);
    };
    this.id = this.board.addChart(this);
    return this.elements;
};
JXG.Chart.prototype = new JXG.GeometryElement;

JXG.Chart.prototype.drawLine = function(board, parents, attributes) {
    var _x = parents[0],
        _y = parents[1];

    /*
    var Fy = function(x) {
        var i, j = 0;

        for(i=0; i<_x.length; i++) {
            if(_x[i] == x) {
                j = i;
                break;
            }
        }

        if(typeof _y[j] == 'function')
            return _y[j]();
        else
            return _y[j];
    }
    */
    
    // not needed
    attributes['fillColor'] = 'none';
    attributes['highlightFillColor'] = 'none';

    var c = board.create('curve', [_x, _y], attributes);
    this.rendNode = c.rendNode;  // This is needed in setProperty
    return c;
};

JXG.Chart.prototype.drawSpline = function(board, parents, attributes) {
    var x = parents[0],
        y = parents[1],
        i;

    // not needed
    attributes['fillColor'] = 'none';
    attributes['highlightFillColor'] = 'none';

    var c = board.create('spline', [x, y], attributes);
    this.rendNode = c.rendNode;  // This is needed in setProperty
    return c;
};

JXG.Chart.prototype.drawFit = function(board, parents, attributes) {
    var x = parents[0],
        y = parents[1],
        deg = (((typeof attributes.degree == 'undefined') || (parseInt(attributes.degree) == NaN)|| (parseInt(attributes.degree) < 1)) ? 1 : parseInt(attributes.degree));

    // not needed
    attributes['fillColor'] = 'none';
    attributes['highlightFillColor'] = 'none';

    var regression = JXG.Math.Numerics.regressionPolynomial(deg, x, y);
    var c = board.create('functiongraph', [regression], attributes);
    this.rendNode = c.rendNode;  // This is needed in setProperty
    return c;
};

JXG.Chart.prototype.drawBar = function(board, parents, attributes) {
    var i, pols = [], x = parents[0], y = parents[1], w, xp0,xp1,xp2, yp, ypL, colorArray, p = [], fill;
    if (attributes['fillOpacity'] == undefined) {
        attributes['fillOpacity'] = 0.6;
    }
    
    // Determine the width of the bars
    if (attributes && attributes['width']) {  // width given
        w = attributes['width'];
    } else {
        if (x.length<=1) {
            w = 1;
        } else {
            // Find minimum distance between to bars.
            w = x[1]-x[0];
            for (i=1;i<x.length-1;i++) {  
                w = (x[i+1]-x[i]<w)?(x[i+1]-x[i]):w;
            }
        }
        w *=0.8;
    }

    fill = attributes['fillColor']
    for (i=0;i<x.length;i++) {        
        if (JXG.isFunction(x[i])) {  // Not yet
            xp0 = function() { return x[i]()-w*0.5; };
            xp1 = function() { return x[i](); };
            xp2 = function() { return x[i]()+w*0.5; };
        } else {
            xp0 = x[i]-w*0.5;
            xp1 = x[i];
            xp2 = x[i]+w*0.5;
        }
        if (JXG.isFunction(y[i])) {  // Not yet
            ypL = yp; //function() { return y[i]()*1.1; };
        } else {
            ypL = y[i]+0.2;
        }
        yp = y[i];
       
        if (attributes['dir']=='horizontal') {  // horizontal bars
            p[0] = board.create('point',[0,xp0], {name:'',fixed:true,visible:false});
            p[1] = board.create('point',[yp,xp0], {name:'',fixed:true,visible:false});
            p[2] = board.create('point',[yp,xp2], {name:'',fixed:true,visible:false});
            p[3] = board.create('point',[0,xp2], {name:'',fixed:true,visible:false});
            if (attributes['labels'] && attributes['labels'][i]) {
                board.create('text',[yp,xp2,attributes['labels'][i]],attributes);
            }
        } else { // vertical bars
            p[0] = board.create('point',[xp0,0], {name:'',fixed:true,visible:false});
            p[1] = board.create('point',[xp0,yp], {name:'',fixed:true,visible:false});
            p[2] = board.create('point',[xp2,yp], {name:'',fixed:true,visible:false});
            p[3] = board.create('point',[xp2,0], {name:'',fixed:true,visible:false});
            if (attributes['labels'] && attributes['labels'][i]) {
                board.create('text',[xp2,yp,attributes['labels'][i]],attributes);
            }
        }
        attributes['withLines'] = false;

        if(typeof fill == 'undefined' && fill == null) {
            colorArray = attributes['colorArray'] || ['#B02B2C','#3F4C6B','#C79810','#D15600','#FFFF88','#C3D9FF','#4096EE','#008C00'];
            attributes['fillColor'] = colorArray[i%colorArray.length];
        }
        pols[i] = board.create('polygon',p,attributes);
    }
    this.rendNode = pols[0].rendNode;  // This is needed in setProperty

    return pols; //[0];  // Not enough! We need pols, but this gives an error in board.setProperty.
};

JXG.Chart.prototype.drawPoints = function(board, parents, attributes) {
    var i;
    var points = [];
    attributes['fixed'] = true;
    attributes['name'] = '';
    var x = parents[0];
    var y = parents[1];
    
    for (i=0;i<x.length;i++) {
        points[i] = board.create('point',[x[i],y[i]], attributes);
    }
    this.rendNode = points[0].rendNode;
    return points; //[0];  // Not enough! We need points, but this gives an error in board.setProperty.
};

JXG.Chart.prototype.drawPie = function(board, parents, attributes) {  // Only 1 array possible as argument 
    var y = parents[0];
    if (y.length<=0) { return; }
    if (typeof y[0] == 'function') { return; } // functions not yet possible

    
    var i;
    var p = [];
    var line = [];
    var arc = [];
    var s = JXG.Math.Statistics.sum(y);
    var colorArray = attributes['colorArray'] || ['#B02B2C','#3F4C6B','#C79810','#D15600','#FFFF88','#C3D9FF','#4096EE','#008C00'];
    var highlightColorArray = attributes['highlightColorArray'] || ['#FF7400'];
    var la = new Array(y.length);
    for(i=0; i<y.length; i++) {
        la[i] = '';
    }
    var labelArray = attributes['labelArray'] || la;
    var radius = attributes['radius'] || 4;
    var myAtts = {};
    if (typeof attributes['highlightOnSector']  =='undefined') {
        attributes['highlightOnSector'] = false;
    }    
    myAtts['name'] = attributes['name'];
    myAtts['id'] = attributes['id'];
    myAtts['strokeWidth'] = attributes['strokeWidth'] || 1;
    myAtts['strokeColor'] = attributes['strokeColor'] || 'none';
    myAtts['straightFirst'] = false;
    myAtts['straightLast'] = false;
    myAtts['fillColor'] = attributes['fillColor'] || '#FFFF88';
    myAtts['fillOpacity'] = attributes['fillOpacity'] || 0.6;
    myAtts['highlightFillColor'] = attributes['highlightFillColor'] || '#FF7400';
    myAtts['highlightStrokeColor'] = attributes['highlightStrokeColor'] || '#FFFFFF';
    myAtts['gradient'] = attributes['gradient'] || 'none';
    var cent = attributes['center'] || [0,0];
    var xc = cent[0];
    var yc = cent[1];

    var center = board.create('point',[xc,yc], {name:'',fixed:true,visible:false});
    p[0] = board.create('point',[radius+xc,0+yc], {name:'',fixed:true,visible:false});
    var rad = 0.0;
    for (i=0;i<y.length;i++) {
        rad += (s!=0)?(2*Math.PI*y[i]/s):0;
        var xcoord = radius*Math.cos(rad)+xc;
        var ycoord = radius*Math.sin(rad)+yc;
        p[i+1] = board.create('point',[xcoord,ycoord], {name:'',fixed:true,visible:false,withLabel:false});
        line[i] = board.create('line',[center,p[i]], 
            {strokeColor:myAtts['strokeColor'], straightFirst:false, straightLast:false, strokeWidth:myAtts['strokeWidth'], strokeOpacity:1.0,withLabel:false,highlightStrokeColor:myAtts['highlightStrokeColor']});
        myAtts['fillColor'] = colorArray[i%colorArray.length];
        myAtts['name'] = labelArray[i];
        if(myAtts['name'] != '') {
            myAtts['withLabel'] = true;
        }
        else {
            myAtts['withLabel'] = false;
        }
        myAtts['labelColor'] = colorArray[i%colorArray.length];
        myAtts['highlightfillColor'] = highlightColorArray[i%highlightColorArray.length];
        arc[i] = board.create('arc',[center,p[i],p[i+1]], myAtts);
        
        if(attributes['highlightOnSector']) {
            arc[i].hasPoint = arc[i].hasPointSector; // overwrite hasPoint so that the whole sector is used for highlighting
        }

    }
    for (i=0;i<y.length;i++) {    
        arc[i].additionalLines = [line[i],line[(i+1)%y.length]];
    }
    this.rendNode = arc[0].rendNode;
    return {arcs:arc, lines:line, points:p, midpoint:center}; //[0];  // Not enough! We need points, but this gives an error in board.setProperty.
};

/**
 * Then, the update function of the renderer
 * is called.  Since a chart is only an abstract element,
 * containing other elements, this function is empty.
 */
JXG.Chart.prototype.updateRenderer = function () {};

/**
 * Update of the defining points
 */
JXG.Chart.prototype.update = function () {
    if (this.needsUpdate) {
        this.updateDataArray();
    }
};

/**
  * For dynamic charts update
  * can be used to compute new entries
  * for the arrays this.dataX and
  * this.dataY. It is used in @see update.
  * Default is an empty method, can be overwritten
  * by the user.
  */
JXG.Chart.prototype.updateDataArray = function () {};


JXG.createChart = function(board, parents, attributes) {
    if((parents.length == 1) && (typeof parents[0] == 'string')) {
        var table = document.getElementById(parents[0]),
            data, row, i, j, col, cell, charts = [], w, x, showRows,
            originalWidth, name, strokeColor, fillColor, hStrokeColor, hFillColor, len;
        if(typeof table != 'undefined') {
            // extract the data
            attributes = JXG.checkAttributes(attributes,{withHeader:true});
            
            table = (new JXG.DataSource()).loadFromTable(parents[0], attributes['withHeader'], attributes['withHeader']);
            data = table.data;
            col = table.columnHeader;
            row = table.rowHeader;

            originalWidth = attributes['width'];
            name = attributes['name'];
            strokeColor = attributes['strokeColor'];
            fillColor = attributes['fillColor'];
            hStrokeColor = attributes['highlightStrokeColor'];
            hFillColor = attributes['highlightFillColor'];

            board.suspendUpdate();

            len = data.length;
            showRows = [];
            if (attributes['rows'] && JXG.isArray(attributes['rows'])) {
                for(i=0; i<len; i++) {
                    for(j=0; j<attributes['rows'].length; j++) {
                        if((attributes['rows'][j] == i) || (attributes['withHeaders'] && attributes['rows'][j] == row[i])) {
                            showRows.push(data[i]);
                            break;
                        }
                    }
                }
            } else {
                showRows = data;
            }

            len = showRows.length;

            for(i=0; i<len; i++) {

                x = [];
                if(attributes['chartStyle'] && attributes['chartStyle'].indexOf('bar') != -1) {
                    if(originalWidth) {
                        w = originalWidth;
                    } else {
                        w = 0.8;
                    }
                    x.push(1 - w/2. + (i+0.5)*w/(1.0*len));
                    for(j=1; j<showRows[i].length; j++) {
                        x.push(x[j-1] + 1);
                    }
                    attributes['width'] = w/(1.0*len);
                }
                
                if(name && name.length == len)
                    attributes['name'] = name[i];
                else if(attributes['withHeaders'])
                    attributes['name'] = col[i];
                
                if(strokeColor && strokeColor.length == len)
                    attributes['strokeColor'] = strokeColor[i];
                else
                    attributes['strokeColor'] = JXG.hsv2rgb(((i+1)/(1.0*len))*360,0.9,0.6);
                
                if(fillColor && fillColor.length == len)
                    attributes['fillColor'] = fillColor[i];
                else
                    attributes['fillColor'] = JXG.hsv2rgb(((i+1)/(1.0*len))*360,0.9,1.0);
                
                if(hStrokeColor && hStrokeColor.length == len)
                    attributes['highlightStrokeColor'] = hStrokeColor[i];
                else
                    attributes['highlightStrokeColor'] = JXG.hsv2rgb(((i+1)/(1.0*len))*360,0.9,1.0);
                
                if(hFillColor && hFillColor.length == len)
                    attributes['highlightFillColor'] = hFillColor[i];
                else
                    attributes['highlightFillColor'] = JXG.hsv2rgb(((i+1)/(1.0*len))*360,0.9,0.6);
                
                if(attributes['chartStyle'] && attributes['chartStyle'].indexOf('bar') != -1) {
                    charts.push(new JXG.Chart(board, [x, showRows[i]], attributes));
                } else
                    charts.push(new JXG.Chart(board, [showRows[i]], attributes));
            }

            board.unsuspendUpdate();

        }
        return charts;
    } else     
        return new JXG.Chart(board, parents, attributes);
};    

JXG.JSXGraph.registerElement('chart', JXG.createChart);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview This file contains code for transformations of geometrical objects. 
 * @author graphjs
 * @version 0.1
 *
 * Possible types:
 * - translate
 * - scale
 * - reflect
 * - rotate
 * - shear
 * - generic
 *
 * Rotation matrix:
 * ( 1    0           0   )
 * ( 0    cos(a)   -sin(a))
 * ( 0    sin(a)   cos(a) )
 *
 * Translation matrix:
 * ( 1  0  0)
 * ( a  1  0)
 * ( b  0  1)

 */
JXG.Transformation = function(board,type, params) { 
    this.elementClass = JXG.OBJECT_CLASS_OTHER;                
    this.matrix = [[1,0,0],[0,1,0],[0,0,1]];
    this.board = board;
    this.isNumericMatrix = false;
    this.setMatrix(board,type,params);
};
JXG.Transformation.prototype = {};

JXG.Transformation.prototype.update = function(){};

/**
 * Set the transformation matrix for different 
 * types of standard transforms
 */
JXG.Transformation.prototype.setMatrix = function(board,type,params) {
    var i;
    
    this.isNumericMatrix = true;
    for (i=0;i<params.length;i++) {
        if (typeof params[i]!='number') {
            this.isNumericMatrix = false;
            break;
        }
    }
    
    if (type=='translate') {
        this.evalParam = JXG.createEvalFunction(board,params,2);
        this.update = function() {
            this.matrix[1][0] = this.evalParam(0);
            this.matrix[2][0] = this.evalParam(1);
        };
    } else if (type=='scale') {
        this.evalParam = JXG.createEvalFunction(board,params,2);
        this.update = function() {
            this.matrix[1][1] = this.evalParam(0); // x
            this.matrix[2][2] = this.evalParam(1); // y
        };
    } else if (type=='reflect') {  // Input: line or two points
        if (params.length<4) { // line or two points
            params[0] = JXG.getReference(board,params[0]);
        }
        if (params.length==2) { // two points
            params[1] = JXG.getReference(board,params[1]);
        }
        if (params.length==4) { // 4 coordinates [px,py,qx,qy]
            this.evalParam = JXG.createEvalFunction(board,params,4);
        }
        this.update = function() {
            var x, y, xoff, yoff, d;
            
            if (params.length==1) { // line
                x = params[0].point2.X()-params[0].point1.X();
                y = params[0].point2.Y()-params[0].point1.Y();
                xoff = params[0].point1.X();
                yoff = params[0].point1.Y();
            } else if (params.length==2){ // two points
                x = params[1].X()-params[0].X();
                y = params[1].Y()-params[0].Y();
                xoff = params[0].X();
                yoff = params[0].Y();
            } else if (params.length==4){ // two points coordinates [px,py,qx,qy]
                x = this.evalParam(2)-this.evalParam(0);
                y = this.evalParam(3)-this.evalParam(1);
                xoff = this.evalParam(0);
                yoff = this.evalParam(1);
            }
            d = x*x+y*y;
            this.matrix[1][1] = (x*x-y*y)/d;
            this.matrix[1][2] = 2*x*y/d;
            this.matrix[2][1] = 2*x*y/d;
            this.matrix[2][2] = (-x*x+y*y)/d;
            this.matrix[1][0] = xoff*(1-this.matrix[1][1])-yoff*this.matrix[1][2];
            this.matrix[2][0] = yoff*(1-this.matrix[2][2])-xoff*this.matrix[2][1];
        };
    } else if (type=='rotate') {
        if (params.length==3) { // angle, x, y
            this.evalParam = JXG.createEvalFunction(board,params,3);
        } else if (params.length<=2) { // angle, p or angle
            this.evalParam = JXG.createEvalFunction(board,params,1);
            if (params.length==2) {
                params[1] = JXG.getReference(board,params[1]);
            } 
        }
        this.update = function() {
            var beta = this.evalParam(0), x, y;
            this.matrix[1][1] = Math.cos(beta); 
            this.matrix[1][2] = -Math.sin(beta);  
            this.matrix[2][1] = Math.sin(beta); 
            this.matrix[2][2] = Math.cos(beta); 
            if (params.length>1) {  // rotate around [x,y] otherwise rotate around [0,0]
                if (params.length==3) {
                    x = this.evalParam(1);
                    y = this.evalParam(2);
                } else {
                    x = params[1].X();
                    y = params[1].Y();
                }
                this.matrix[1][0] = x*(1-Math.cos(beta))+y*Math.sin(beta);
                this.matrix[2][0] = y*(1-Math.cos(beta))-x*Math.sin(beta);
            }
        };
    } else if (type=='shear') {
        this.evalParam = JXG.createEvalFunction(board,params,1);
        this.update = function() {
            var beta = this.evalParam(0);
            this.matrix[1][1] = Math.tan(beta); 
        };
    } else if (type=='generic') {
        this.evalParam = JXG.createEvalFunction(board,params,9);
        this.update = function() {
            this.matrix[0][0] = this.evalParam(0); 
            this.matrix[0][1] = this.evalParam(1); 
            this.matrix[0][2] = this.evalParam(2); 
            this.matrix[1][0] = this.evalParam(3); 
            this.matrix[1][1] = this.evalParam(4); 
            this.matrix[1][2] = this.evalParam(5); 
            this.matrix[2][0] = this.evalParam(6); 
            this.matrix[2][1] = this.evalParam(7); 
            this.matrix[2][2] = this.evalParam(8); 
        };
    }
};

/**
 * Transform a GeometryElement:
 * First, update the matrix
 * Second, do the matrix-vector-multiplication
 *
 * @param {JXG.GeometryElement} element, which is transformed
 */
JXG.Transformation.prototype.apply = function(p){
    this.update();
    if (arguments[1]!=null) {
        return JXG.Math.matVecMult(this.matrix,p.initialCoords.usrCoords);
    } else {
        return JXG.Math.matVecMult(this.matrix,p.coords.usrCoords);
    }
};

/**
 * Apply a transformation once to a GeometryElement.
 * If it is a free point, then it can be dragged around later
 * and will overwrite the transformed coordinates.
 */
JXG.Transformation.prototype.applyOnce = function(p){
    var c, len, i;
    if (!JXG.isArray(p)) {   
        this.update();
        c = JXG.Math.matVecMult(this.matrix,p.coords.usrCoords);
        p.coords.setCoordinates(JXG.COORDS_BY_USER,[c[1],c[2]]);
    } else {
        len = p.length;
        for (i=0; i<len; i++) {
            this.update();
            c = JXG.Math.matVecMult(this.matrix,p[i].coords.usrCoords);
            p[i].coords.setCoordinates(JXG.COORDS_BY_USER,[c[1],c[2]]);
        }
    }
};

/**
 * Bind a transformation to a GeometryElement
 */
JXG.Transformation.prototype.bindTo = function(p){
    var i, len;
    if (JXG.isArray(p)) {   
        len = p.length;
        for (i=0; i<len; i++) {
            p[i].transformations.push(this);
        }
    } else {
        p.transformations.push(this);
    }
};

JXG.Transformation.prototype.setProperty = function(term) {};

/**
 * Multiplication of a transformation t from the right.
 * this = t join this
 */
JXG.Transformation.prototype.melt = function(t){
    var res = [], i, len, len0, k, s, j;
    
    len = t.matrix.length;
    len0 = this.matrix[0].length;
    
    for (i=0;i<len;i++) {
        res[i] = [];
    }
    this.update();
    t.update();
    for (i=0;i<len;i++) {
        for (j=0;j<len0;j++) {
            s = 0;
            for (k=0;k<len;k++) {
                s += t.matrix[i][k]*this.matrix[k][j];
            }
            res[i][j] = s;
        }
    }
    this.update = function() {
        var len = this.matrix.length,
            len0 = this.matrix[0].length;
        for (i=0;i<len;i++) {
            for (j=0;j<len0;j++) {
                this.matrix[i][j] = res[i][j];
            }
        }
    };
    return true;
};

JXG.createTransform = function(board, parentArr, atts) {
    return new JXG.Transformation(board,atts['type'],parentArr);
};

JXG.JSXGraph.registerElement('transform', JXG.createTransform);
 /*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * @fileoverview The JSXGraph object Turtle is defined. It acts like
 * "turtle graphics".
 * @author A.W.
 */

/**
 * Constructs a new Turtle object.
 * @class This is the Turtle class. 
 * It is derived from {@link JXG.GeometryElement}.
 * It stores all properties required
 * to move a turtle.
 * @constructor
 * @param {String} JXG.board The board the new turtle is drawn on.
 * @param {Array}  [x,y,angle] Start position and start direction of the turtle. Possible values are
 * [x,y,angle]
 * [[x,y],angle]
 * [x,y]
 * [[x,y]]
 * @param {Object} attributes Attributes to change the visual properties of the turtle object
 * All angles are in degrees.
  */
JXG.Turtle = function (board, parents, attributes) {
    var x,y,dir;
    this.type = JXG.OBJECT_TYPE_TURTLE;
    this.turtleIsHidden = false;
    this.board = board;
    this.attributes = JXG.checkAttributes(attributes,{withLabel:false,layer:null});
    this.attributes.straightFirst = false;
    this.attributes.straightLast = false;
    x = 0;
    y = 0;
    dir = 90;
    if (parents.length!=0) {
        if (parents.length==3) {   // [x,y,dir]
            // Only numbers are accepted at the moment
            x = parents[0];
            y = parents[1];
            dir = parents[2];
        } else if (parents.length==2) {
            if (JXG.isArray(parents[0])) {  // [[x,y],dir]
                x = parents[0][0];
                y = parents[0][1];
                dir = parents[1];
            } else {  // [x,y]
                x = parents[0];
                y = parents[1];
            }
        } else { // [[x,y]]
            x = parents[0][0];
            y = parents[0][1];
        }
    }
    
    this.init(x,y,dir);
    return this;
};
JXG.Turtle.prototype = new JXG.GeometryElement;

/**
* Initialize a new turtle or reinitialize  a turtle after {@link #clearscreen}.
* @private
*/
JXG.Turtle.prototype.init = function(x,y,dir) {
    this.arrowLen = 20.0/Math.sqrt(this.board.unitX*this.board.unitX+this.board.unitY*this.board.unitY);

    this.pos = [x,y];
    this.isPenDown = true;
    this.dir = 90;
    this.stack = [];
    this.objects = [];
    this.attributes.curveType = 'plot';
    this.curve = this.board.create('curve',[[this.pos[0]],[this.pos[1]]],this.attributes);
    this.objects.push(this.curve);

    this.turtle = this.board.create('point',this.pos,{fixed:true,name:' ',visible:false,withLabel:false});
    this.objects.push(this.turtle);
    
    this.turtle2 = this.board.create('point',[this.pos[0],this.pos[1]+this.arrowLen],
            {fixed:true,name:' ',visible:false,withLabel:false});
    this.objects.push(this.turtle2);
    
    var w = this.attributes.strokeWidth || this.attributes.strokewidth || 2;  // Attention; should be moved to Options.js
    this.arrow = this.board.create('line',[this.turtle,this.turtle2],
            {lastArrow:true,strokeColor:'#ff0000',straightFirst:false,straightLast:false,strokeWidth:w,withLabel:false});
    this.objects.push(this.arrow);

    this.right(90-dir);
    this.board.update();
};

/**
* Move the turtle forward.
* @param {float} length of forward move in user coordinates
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.forward = function(len) {
    if (len==0) { return; }
    var dx = len*Math.cos(this.dir*Math.PI/180.0);
    var dy = len*Math.sin(this.dir*Math.PI/180.0);
    if (!this.turtleIsHidden) {
        var t = this.board.create('transform', [dx,dy], {type:'translate'});
        t.applyOnce(this.turtle);
        t.applyOnce(this.turtle2);
    }
    if (this.isPenDown) if (this.curve.dataX.length>=8192) { // IE workaround
        this.curve = this.board.create('curve',
               [[this.pos[0]],[this.pos[1]]],this.attributes);
        this.objects.push(this.curve);
    }
    this.pos[0] += dx;
    this.pos[1] += dy;
    if (this.isPenDown) {
        this.curve.dataX.push(this.pos[0]);
        this.curve.dataY.push(this.pos[1]);
    }
    this.board.update();
    return this;
};
     
/**
* Move the turtle backwards.
* @param {float} length of backwards move in user coordinates
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.back = function(len) {
    return this.forward(-len);
};
     
/**
* Rotate the turtle direction to the right
* @param {float} angle of the rotation in degrees
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.right = function(angle) {
    this.dir -= angle;
    this.dir %= 360.0;
    if (!this.turtleIsHidden) {
        var t = this.board.create('transform', [-angle*Math.PI/180.0,this.turtle], {type:'rotate'});
        t.applyOnce(this.turtle2);
    }
    this.board.update();
    return this;
};
     
/**
* Rotate the turtle direction to the right.
* @param {float} angle of the rotation in degrees
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.left = function(angle) {
    return this.right(-angle);
};

/**
* Pen up, stops visible drawing
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.penUp = function() {
    this.isPenDown = false;
    return this;
};

/**
* Pen down, continues visible drawing
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.penDown = function() {
    this.isPenDown = true;
    this.curve = this.board.create('curve',[[this.pos[0]],[this.pos[1]]],this.attributes);
    this.objects.push(this.curve);
            
    return this;
};

/**
*  Removes the turtle curve from the board. The turtle stays in its position.
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.clean = function() {
    for(var i=0;i<this.objects.length;i++) {
        var el = this.objects[i];
        if (el.type==JXG.OBJECT_TYPE_CURVE) {
            this.board.removeObject(el.id);
            this.objects.splice(i,1);
        }
    }
    this.curve = this.board.create('curve',
              [[this.pos[0]],[this.pos[1]]],this.attributes);
    this.objects.push(this.curve);
    this.board.update();
    return this;
};

/**
*  Removes the turtle completely and resets it to its initial position and direction.
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.clearScreen = function() {
    for(var i=0;i<this.objects.length;i++) {
        var el = this.objects[i];
        this.board.removeObject(el.id);
    }
    this.init(0,0,90);
    return this;
};

/**
*  Moves the turtle without drawing to a new position
* @param {float} x new x- coordinate 
* @param {float} y new y- coordinate 
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.setPos = function(x,y) {
    if (JXG.isArray(x)) {
        this.pos = x;
    } else {
        this.pos = [x,y];
    }
    if (!this.turtleIsHidden) {
        this.turtle.setPositionDirectly(JXG.COORDS_BY_USER,x,y);
        this.turtle2.setPositionDirectly(JXG.COORDS_BY_USER,x,y+this.arrowLen);
        var t = this.board.create('transform', 
                [-(this.dir-90)*Math.PI/180.0,this.turtle], {type:'rotate'});
        t.applyOnce(this.turtle2);
    }
    this.curve = this.board.create('curve',[[this.pos[0]],[this.pos[1]]],this.attributes);
    this.objects.push(this.curve);
    this.board.update();
    return this;
};

/**
*  Sets the pen size. Equivalent to setProperty({strokeWidth:size})
* @param {float} size
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.setPenSize = function(size) { 
    this.attributes.strokeWidth = size; 
    this.curve = this.board.create('curve',[[this.pos[0]],[this.pos[1]]],this.attributes);
    this.objects.push(this.curve);
    return this;
};

/**
*  Sets the pen color. Equivalent to setProperty({strokeColor:color})
* @param {string} color
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.setPenColor = function(colStr) { 
    this.attributes.strokeColor = colStr; 
    this.curve = this.board.create('curve',[[this.pos[0]],[this.pos[1]]],this.attributes);
    this.objects.push(this.curve);
    return this;
};

/**
*  Sets the highlight pen color. Equivalent to setProperty({highlightStrokeColor:color})
* @param {string} color
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.setHighlightPenColor = function(colStr) { 
    this.attributes.highlightStrokeColor = colStr; 
    this.curve = this.board.create('curve',[[this.pos[0]],[this.pos[1]]],this.attributes);
    this.objects.push(this.curve);
    return this;
};

/**
* Sets properties of the turtle, see also {@link JXG.GeometryElement#setProperty}.
* Sets the property for all curves of the turtle.
* @param {Object} key:value pairs
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.setProperty = function() {
    var pair;
    var pairRaw;
    var i, el;
    var key;
    for (i=0; i<arguments.length; i++) {
        pairRaw = arguments[i];
        if (typeof pairRaw == 'string') {    // pairRaw is string of the form 'key:value'
            pair = pairRaw.split(':');
        } else if (!JXG.isArray(pairRaw)) {    
            // pairRaw consists of objects of the form {key1:value1,key2:value2,...}
            for (var key in pairRaw) {
                this.setProperty([key,pairRaw[key]]);
            }
            return this;
        } else {                             // pairRaw consists of array [key,value]
            pair = pairRaw;
        }
        this.attributes[pair[0]] = pair[1];
    }
    for (i=0; i<this.objects.length; i++) {
        el = this.objects[i];
        if (el.type==JXG.OBJECT_TYPE_CURVE) {
            el.setProperty(this.attributes);
        }
    }
    //this.curve = this.board.create('curve',[[this.pos[0]],[this.pos[1]]],this.attributes);
    //this.objects.push(this.curve);
    return this;
};

/**
*  Sets the visibility of the turtle head to true,
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.showTurtle = function() { 
    this.turtleIsHidden = false; 
    this.arrow.setProperty('visible:true');
    this.setPos(this.pos[0],this.pos[1]);
    this.board.update();
    return this;
};

/**
*  Sets the visibility of the turtle head to false,
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.hideTurtle = function() { 
    this.turtleIsHidden = true;
    this.arrow.setProperty('visible:false');
    this.setPos(this.pos[0],this.pos[1]);
    this.board.update();
    return this;
};

/**
*  Moves the turtle to position [0,0].
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.home = function() { 
    this.pos = [0,0];
    this.setPos(this.pos[0],this.pos[1]);
    return this;
};

/**
*  Pushes the position of the turtle on the stack.
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.pushTurtle = function() { 
    this.stack.push([this.pos[0],this.pos[1],this.dir]);
    return this;
};

/**
*  Gets the last position of the turtle on the stack, sets the turtle to this position and removes this 
* position from the stack.
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.popTurtle = function() { 
    var status = this.stack.pop();
    this.pos[0] = status[0];
    this.pos[1] = status[1];
    this.dir = status[2];
    this.setPos(this.pos[0],this.pos[1]);
    return this;
};

/**
* Rotates the turtle into a new direction.
* There are two possibilities:
* @param {float} angle New direction to look to
* or
* @param {float} x New x coordinate to look to
* @param {float} y New y coordinate to look to
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.lookTo = function(target) { 
    if (JXG.isArray(target)) {
        var ax = this.pos[0];
        var ay = this.pos[1];
        var bx = target[0];
        var by = target[1];
        var beta; 
        // Rotate by the slope of the line [this.pos, target]
        var sgn = (bx-ax>0)?1:-1;
        if (Math.abs(bx-ax)>0.0000001) {
            beta = Math.atan2(by-ay,bx-ax)+((sgn<0)?Math.PI:0);  
        } else {
            beta = ((by-ay>0)?0.5:-0.5)*Math.PI;
        }
        this.right(this.dir-(beta*180/Math.PI));
    } else if (JXG.isNumber(target)) {
        this.right(this.dir-(target));
    }
    return this;
};

/**
* Moves the turtle to a given coordinate pair.
* The direction is not changed.
* @param {float} x New x coordinate to look to
* @param {float} y New y coordinate to look to
* @type {JXG.Turtle}
* @return pointer to the turtle object
*/
JXG.Turtle.prototype.moveTo = function(target) { 
    if (JXG.isArray(target)) {
        var dx = target[0]-this.pos[0];
        var dy = target[1]-this.pos[1];
        if (!this.turtleIsHidden) {
            var t = this.board.create('transform', [dx,dy], {type:'translate'});
            t.applyOnce(this.turtle);
            t.applyOnce(this.turtle2);
        }
        if (this.isPenDown) if (this.curve.dataX.length>=8192) { // IE workaround
            this.curve = this.board.create('curve',
               [[this.pos[0]],[this.pos[1]]],this.attributes);
            this.objects.push(this.curve);
        }
        this.pos[0] = target[0];
        this.pos[1] = target[1];
        if (this.isPenDown) {
            this.curve.dataX.push(this.pos[0]);
            this.curve.dataY.push(this.pos[1]);
        }
        this.board.update();
    }
    return this;
};

/**
  * Alias for {@link #forward}
  */
JXG.Turtle.prototype.fd = function(len) { return this.forward(len); };
/**
  * Alias for {@link #back}
  */
JXG.Turtle.prototype.bk = function(len) { return this.back(len); };
/**
  * Alias for {@link #left}
  */
JXG.Turtle.prototype.lt = function(angle) { return this.left(angle); };
/**
  * Alias for {@link #right}
  */
JXG.Turtle.prototype.rt = function(angle) { return this.right(angle); };
/**
  * Alias for {@link #penUp}
  */
JXG.Turtle.prototype.pu = function() { return this.penUp(); };
/**
  * Alias for {@link #penDown}
  */
JXG.Turtle.prototype.pd = function() { return this.penDown(); };
/**
  * Alias for {@link #hideTurtle}
  */
JXG.Turtle.prototype.ht = function() { return this.hideTurtle(); };
/**
  * Alias for {@link #showTurtle}
  */
JXG.Turtle.prototype.st = function() { return this.showTurtle(); };
/**
  * Alias for {@link #clearScreen}
  */
JXG.Turtle.prototype.cs = function() { return this.clearScreen(); };
/**
  * Alias for {@link #pushTurtle}
  */
JXG.Turtle.prototype.push = function() { return this.pushTurtle(); };
/**
  * Alias for {@link #popTurtle}
  */
JXG.Turtle.prototype.pop = function() { return this.popTurtle(); };

/**
* @return x-coordinate of the turtle position
* @type {float}
*/
JXG.Turtle.prototype.X = function(target) { 
    return this.pos[0]; //this.turtle.X();
};

/**
* @return y-coordinate of the turtle position
* @type {float}
*/
JXG.Turtle.prototype.Y = function(target) { 
    return this.pos[1]; //this.turtle.Y();
};

/**
 * Checks whether (x,y) is near the curve.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @param {y} Find closest point on the curve to (x,y)
 * @return {bool} True if (x,y) is near the curve, False otherwise.
 */
JXG.Turtle.prototype.hasPoint = function (x,y) {
    var i, el;
    for(i=0;i<this.objects.length;i++) {  // run through all curves of this turtle
        el = this.objects[i];
        if (el.type==JXG.OBJECT_TYPE_CURVE) {
            if (el.hasPoint(x,y)) {
                return true;              // So what??? All other curves have to be notified now (for highlighting)
                                          // This has to be done, yet.
            }
        }
    }
    return false;
};

/**
 * Creates a new turtle
 * @param {JXG.Board} board The board the turtle is put on.
 * @param {Array} parents 
 * @param {Object} attributs Object containing properties for the element such as stroke-color and visibility. See {@link JXG.GeometryElement#setProperty}
 * @type JXG.Turtle
 * @return Reference to the created turtle object.
 */
JXG.createTurtle = function(board, parents, attributes) {
    if (parents==null) {
        var parents = [];
    }
    return new JXG.Turtle(board,parents,attributes);
};

JXG.JSXGraph.registerElement('turtle', JXG.createTurtle);

/**
 * Functions for color conversions. Based on a class to parse color values by Stoyan Stefanov <sstoo@gmail.com>
 * @see http://www.phpied.com/rgb-color-parser-in-javascript/
 */

/**
 * Converts a valid HTML/CSS color string into a rgb value array. This is the base
 * function for the following wrapper functions which only adjust the output to
 * different flavors like an object, string or hex values.
 * @parameter {string} color_string A valid HTML or CSS styled color value, e.g. #12ab21, #abc, black, or rgb(12, 132, 233) <strong>or</string>
 * @parameter {array} color_array Array containing three color values either from 0.0 to 1.0 or from 0 to 255. They will be interpreted as red, green, and blue values <strong>OR</strong>
 * @parameter {number} r,g,b Three color values r, g, and b like those in the array variant.
 * @type array
 * @return RGB color values as an array [r, g, b] which component's are between 0 and 255.
 */
JXG.rgbParser = function() {

    if(arguments.length == 0)
        return;

    if(arguments.length >= 3) {
        arguments[0] = [arguments[0], arguments[1], arguments[2]];
        arguments.length = 1;
    }

    var color_string = arguments[0];
    if(JXG.isArray(color_string)) {
        var testFloat = false, i;
        for(i=0; i<3; i++)
            testFloat |= /\./.test(arguments[0][i].toString());
        for(i=0; i<3; i++)
            testFloat &= (arguments[0][i] >= 0.0) & (arguments[0][i] <= 1.0);

        if(testFloat)
            return [Math.ceil(arguments[0][0] * 255), Math.ceil(arguments[0][1] * 255), Math.ceil(arguments[0][2] * 255)];
        else {
            arguments[0].length = 3;
            return arguments[0];
        }
    } else if(typeof arguments[0] == 'string') {
        color_string = arguments[0];
    }

    var r, g, b;

    // strip any leading #
    if (color_string.charAt(0) == '#') { // remove # if any
        color_string = color_string.substr(1,6);
    }

    color_string = color_string.replace(/ /g,'');
    color_string = color_string.toLowerCase();

    // before getting into regexps, try simple matches
    // and overwrite the input
    var simple_colors = {
        aliceblue: 'f0f8ff',
        antiquewhite: 'faebd7',
        aqua: '00ffff',
        aquamarine: '7fffd4',
        azure: 'f0ffff',
        beige: 'f5f5dc',
        bisque: 'ffe4c4',
        black: '000000',
        blanchedalmond: 'ffebcd',
        blue: '0000ff',
        blueviolet: '8a2be2',
        brown: 'a52a2a',
        burlywood: 'deb887',
        cadetblue: '5f9ea0',
        chartreuse: '7fff00',
        chocolate: 'd2691e',
        coral: 'ff7f50',
        cornflowerblue: '6495ed',
        cornsilk: 'fff8dc',
        crimson: 'dc143c',
        cyan: '00ffff',
        darkblue: '00008b',
        darkcyan: '008b8b',
        darkgoldenrod: 'b8860b',
        darkgray: 'a9a9a9',
        darkgreen: '006400',
        darkkhaki: 'bdb76b',
        darkmagenta: '8b008b',
        darkolivegreen: '556b2f',
        darkorange: 'ff8c00',
        darkorchid: '9932cc',
        darkred: '8b0000',
        darksalmon: 'e9967a',
        darkseagreen: '8fbc8f',
        darkslateblue: '483d8b',
        darkslategray: '2f4f4f',
        darkturquoise: '00ced1',
        darkviolet: '9400d3',
        deeppink: 'ff1493',
        deepskyblue: '00bfff',
        dimgray: '696969',
        dodgerblue: '1e90ff',
        feldspar: 'd19275',
        firebrick: 'b22222',
        floralwhite: 'fffaf0',
        forestgreen: '228b22',
        fuchsia: 'ff00ff',
        gainsboro: 'dcdcdc',
        ghostwhite: 'f8f8ff',
        gold: 'ffd700',
        goldenrod: 'daa520',
        gray: '808080',
        green: '008000',
        greenyellow: 'adff2f',
        honeydew: 'f0fff0',
        hotpink: 'ff69b4',
        indianred : 'cd5c5c',
        indigo : '4b0082',
        ivory: 'fffff0',
        khaki: 'f0e68c',
        lavender: 'e6e6fa',
        lavenderblush: 'fff0f5',
        lawngreen: '7cfc00',
        lemonchiffon: 'fffacd',
        lightblue: 'add8e6',
        lightcoral: 'f08080',
        lightcyan: 'e0ffff',
        lightgoldenrodyellow: 'fafad2',
        lightgrey: 'd3d3d3',
        lightgreen: '90ee90',
        lightpink: 'ffb6c1',
        lightsalmon: 'ffa07a',
        lightseagreen: '20b2aa',
        lightskyblue: '87cefa',
        lightslateblue: '8470ff',
        lightslategray: '778899',
        lightsteelblue: 'b0c4de',
        lightyellow: 'ffffe0',
        lime: '00ff00',
        limegreen: '32cd32',
        linen: 'faf0e6',
        magenta: 'ff00ff',
        maroon: '800000',
        mediumaquamarine: '66cdaa',
        mediumblue: '0000cd',
        mediumorchid: 'ba55d3',
        mediumpurple: '9370d8',
        mediumseagreen: '3cb371',
        mediumslateblue: '7b68ee',
        mediumspringgreen: '00fa9a',
        mediumturquoise: '48d1cc',
        mediumvioletred: 'c71585',
        midnightblue: '191970',
        mintcream: 'f5fffa',
        mistyrose: 'ffe4e1',
        moccasin: 'ffe4b5',
        navajowhite: 'ffdead',
        navy: '000080',
        oldlace: 'fdf5e6',
        olive: '808000',
        olivedrab: '6b8e23',
        orange: 'ffa500',
        orangered: 'ff4500',
        orchid: 'da70d6',
        palegoldenrod: 'eee8aa',
        palegreen: '98fb98',
        paleturquoise: 'afeeee',
        palevioletred: 'd87093',
        papayawhip: 'ffefd5',
        peachpuff: 'ffdab9',
        peru: 'cd853f',
        pink: 'ffc0cb',
        plum: 'dda0dd',
        powderblue: 'b0e0e6',
        purple: '800080',
        red: 'ff0000',
        rosybrown: 'bc8f8f',
        royalblue: '4169e1',
        saddlebrown: '8b4513',
        salmon: 'fa8072',
        sandybrown: 'f4a460',
        seagreen: '2e8b57',
        seashell: 'fff5ee',
        sienna: 'a0522d',
        silver: 'c0c0c0',
        skyblue: '87ceeb',
        slateblue: '6a5acd',
        slategray: '708090',
        snow: 'fffafa',
        springgreen: '00ff7f',
        steelblue: '4682b4',
        tan: 'd2b48c',
        teal: '008080',
        thistle: 'd8bfd8',
        tomato: 'ff6347',
        turquoise: '40e0d0',
        violet: 'ee82ee',
        violetred: 'd02090',
        wheat: 'f5deb3',
        white: 'ffffff',
        whitesmoke: 'f5f5f5',
        yellow: 'ffff00',
        yellowgreen: '9acd32'
    };
    for (var key in simple_colors) {
        if (color_string == key) {
            color_string = simple_colors[key];
        }
    }
    // end of simple type-in colors

    // array of color definition objects
    var color_defs = [
        {
            re: /^rgb\((\d{1,3}),\s*(\d{1,3}),\s*(\d{1,3})\)$/,
            example: ['rgb(123, 234, 45)', 'rgb(255,234,245)'],
            process: function (bits){
                return [
                    parseInt(bits[1]),
                    parseInt(bits[2]),
                    parseInt(bits[3])
                ];
            }
        },
        {
            re: /^(\w{2})(\w{2})(\w{2})$/,
            example: ['#00ff00', '336699'],
            process: function (bits){
                return [
                    parseInt(bits[1], 16),
                    parseInt(bits[2], 16),
                    parseInt(bits[3], 16)
                ];
            }
        },
        {
            re: /^(\w{1})(\w{1})(\w{1})$/,
            example: ['#fb0', 'f0f'],
            process: function (bits){
                return [
                    parseInt(bits[1] + bits[1], 16),
                    parseInt(bits[2] + bits[2], 16),
                    parseInt(bits[3] + bits[3], 16)
                ];
            }
        }
    ];

    // search through the definitions to find a match
    for (var i = 0; i < color_defs.length; i++) {
        var re = color_defs[i].re;
        var processor = color_defs[i].process;
        var bits = re.exec(color_string);
        if (bits) {
            channels = processor(bits);
            r = channels[0];
            g = channels[1];
            b = channels[2];
        }

    }

    // validate/cleanup values
    r = (r < 0 || isNaN(r)) ? 0 : ((r > 255) ? 255 : r);
    g = (g < 0 || isNaN(g)) ? 0 : ((g > 255) ? 255 : g);
    b = (b < 0 || isNaN(b)) ? 0 : ((b > 255) ? 255 : b);

    return [r, g, b];
};

/**
 * Returns output of JXG.rgbParser as a CSS styled rgb() string.
 */
JXG.rgb2css = function () {
    var r, g, b;
    r = JXG.rgbParser.apply(JXG.rgbParser, arguments);
    g = r[1];
    b = r[2];
    r = r[0];
    return 'rgb(' + r + ', ' + g + ', ' + b + ')';
};

/**
 * Returns array returned by JXG.rgbParser as a HTML rgb string.
 */
JXG.rgb2hex = function () {
    var r, g, b;
    r = JXG.rgbParser.apply(JXG.rgbParser, arguments);
    g = r[1];
    b = r[2];
    r = r[0];
    r = r.toString(16);
    g = g.toString(16);
    b = b.toString(16);
    if (r.length == 1) r = '0' + r;
    if (g.length == 1) g = '0' + g;
    if (b.length == 1) b = '0' + b;
    return '#' + r + g + b;
};

/**
* Converts HSV color to RGB color.
* Based on C Code in "Computer Graphics -- Principles and Practice,"
* Foley et al, 1996, p. 593.
* See also http://www.efg2.com/Lab/Graphics/Colors/HSV.htm  
* @param {float} H value between 0 and 360
* @param {float} S value between 0.0 (shade of gray) to 1.0 (pure color)
* @param {float} V value between 0.0 (black) to 1.0 (white)
* @return {string} RGB color string
*/
JXG.hsv2rgb = function(H,S,V) {
    var R,G,B, f,i,hTemp, p,q,t;
    H = ((H%360.0)+360.0)%360;
    if (S==0) {
        if (isNaN(H) || H < JXG.Math.eps) {
            R = V;
            G = V;
            B = V;
        } else {
            return '#ffffff';
        }
    } else {
        if (H>=360) {
            hTemp = 0.0;
        } else {
            hTemp = H;
        }
        hTemp = hTemp / 60;     // h is now IN [0,6)
        i = Math.floor(hTemp);        // largest integer <= h
        f = hTemp - i;                  // fractional part of h
        p = V * (1.0 - S);
        q = V * (1.0 - (S * f));
        t = V * (1.0 - (S * (1.0 - f)));
        switch (i) {
            case 0: R = V; G = t;  B = p; break;
            case 1: R = q; G = V;  B = p; break;
            case 2: R = p; G = V;  B = t; break;
            case 3: R = p; G = q;  B = V; break;
            case 4: R = t; G = p;  B = V; break;
            case 5: R = V; G = p;  B = q; break;
        }
    }
    R = Math.round(R*255).toString(16); R = (R.length==2)?R:((R.length==1)?'0'+R:'00');
    G = Math.round(G*255).toString(16); G = (G.length==2)?G:((G.length==1)?'0'+G:'00');
    B = Math.round(B*255).toString(16); B = (B.length==2)?B:((B.length==1)?'0'+B:'00');
    return ['#',R,G,B].join(''); 
};

/**
 * Converts r, g, b color to h, s, v.
 * See http://zach.in.tu-clausthal.de/teaching/cg1_0708/folien/13_color_3_4up.pdf for more information.
 * @param {number} r Amount of red in color. Number between 0 and 255.
 * @param {number} g Amount of green. Number between 0 and 255.
 * @param {number} b Amount of blue. Number between 0 and 255.
 * @type Object
 * @return Hashmap containing h,s, and v field.
 */
JXG.rgb2hsv = function() {
    var r, g, b, fr, fg, fb, fmax, fmin, h, s, v, max, min, stx;
    r = JXG.rgbParser.apply(JXG.rgbParser, arguments);
    g = r[1];
    b = r[2];
    r = r[0];
    stx = JXG.Math.Statistics;
    fr = r/255.;
    fg = g/255.;
    fb = b/255.;
    max = stx.max([r, g, b]);
    min = stx.min([r, g, b]);
    fmax = max/255.;
    fmin = min/255.;

    v = fmax;

    s = 0.;
    if(v>0) {
        s = (v-fmin)/(v*1.);
    }

    h = 1./(fmax-fmin);
    if(s > 0) {
        if(max==r)
            h = (fg-fb)*h;
        else if(max==g)
            h = 2 + (fb-fr)*h;
        else
            h = 4 + (fr-fg)*h;
    }

    h *= 60;
    if(h < 0)
        h += 360;

    if(max==min)
        h = 0.;

    return [h, s, v];
};


/**
 * Convert RGB color information to LMS color space.
 * @param {number} r Amount of red in color. Number between 0 and 255.
 * @param {number} g Amount of green. Number between 0 and 255.
 * @param {number} b Amount of blue. Number between 0 and 255.
 * @type Object
 * @return Hashmap containing the L, M, S cone values.
 */
JXG.rgb2LMS = function() {
    var r, g, b, l, m, s, ret
        // constants
        matrix = [[0.05059983, 0.08585369, 0.00952420], [0.01893033, 0.08925308, 0.01370054], [0.00292202, 0.00975732, 0.07145979]];

    r = JXG.rgbParser.apply(JXG.rgbParser, arguments);
    g = r[1];
    b = r[2];
    r = r[0];

    // de-gamma
    // Maybe this can be made faster by using a cache
    r = Math.pow(r, 0.476190476);
    g = Math.pow(g, 0.476190476);
    b = Math.pow(b, 0.476190476);

    l = r * matrix[0][0] + g * matrix[0][1] + b * matrix[0][2];
    m = r * matrix[1][0] + g * matrix[1][1] + b * matrix[1][2];
    s = r * matrix[2][0] + g * matrix[2][1] + b * matrix[2][2];

    ret = [l, m, s];
    ret.l = l;
    ret.m = m;
    ret.s = s;

    return ret;
};
/**
 * Convert color information from LMS to RGB color space.
 * @param {number} l Amount of l value.
 * @param {number} m Amount of m value.
 * @param {number} s Amount of s value.
 * @type Object
 * @return Hashmap containing the r, g, b values.
 */
JXG.LMS2rgb = function(l, m, s) {
    var r, g, b, ret
        // constants
        matrix = [[30.830854, -29.832659, 1.610474], [-6.481468, 17.715578, -2.532642], [-0.375690, -1.199062, 14.273846]];

    // transform back to rgb
    r = l * matrix[0][0] + m * matrix[0][1] + s * matrix[0][2];
    g = l * matrix[1][0] + m * matrix[1][1] + s * matrix[1][2];
    b = l * matrix[2][0] + m * matrix[2][1] + s * matrix[2][2];

    // re-gamma, inspired by GIMP modules/display-filter-color-blind.c:
    // Copyright (C) 2002-2003 Michael Natterer <mitch@gimp.org>,
    //                         Sven Neumann <sven@gimp.org>,
    //                         Robert Dougherty <bob@vischeck.com> and
    //                         Alex Wade <alex@vischeck.com>
    // This code is an implementation of an algorithm described by Hans Brettel,
    // Francoise Vienot and John Mollon in the Journal of the Optical Society of
    // America V14(10), pg 2647. (See http://vischeck.com/ for more info.)
    lut_lookup = function (value) {
        var offset = 127, step = 64;

        while (step > 0) {
            if (Math.pow(offset, 0.476190476) > value) {
                offset -= step;
            } else {
                if (Math.pow(offset+1, 0.476190476) > value)
                    return offset;

                offset += step;
            }

            step /= 2;
        }

        /*  the algorithm above can't reach 255  */
        if (offset == 254 && 13.994955247 < value)
            return 255;

        return offset;
    };


    r = lut_lookup(r);
    g = lut_lookup(g);
    b = lut_lookup(b);

    ret = [r, g, b];
    ret.r = r;
    ret.g = g;
    ret.b = b;

    return ret;
};

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

JXG.Board.prototype.angle = function(A, B, C){ return this.algebra.angle(A,B,C); };
JXG.Board.prototype.rad = function(A, B, C){ return this.algebra.rad(A,B,C); };
JXG.Board.prototype.distance = function(arr1, arr2){ return this.algebra.distance(arr1,arr2); };
JXG.Board.prototype.pow = function(a, b){ return this.algebra.pow(a,b); };
JXG.Board.prototype.round = function(x, n){ return (x).toFixed(n); };
JXG.Board.prototype.cosh = function(x){ return JXG.Math.cosh(x); };
JXG.Board.prototype.sinh = function(x){ return JXG.Math.sinh(x); };
JXG.Board.prototype.sgn = function(x) { return (x==0 ? 0 : x/(Math.abs(x))); };
JXG.Board.prototype.D = function(f,obj){ return JXG.Math.Numerics.D(f,obj); };
JXG.Board.prototype.I = function(interval,f){ return JXG.Math.Numerics.I(interval,f); };
JXG.Board.prototype.root = function(f,x,obj){ return JXG.Math.Numerics.root(f,x,obj); };
JXG.Board.prototype.lagrangePolynomial = function(p){ return JXG.Math.Numerics.lagrangePolynomial(p); };
JXG.Board.prototype.neville = function(p){ return JXG.Math.Numerics.neville(p); };
JXG.Board.prototype.riemannsum = function(f,n,type,start,end){ return JXG.Math.Numerics.riemannsum(f,n,type,start,end); };

JXG.Board.prototype.abs = Math.abs;
JXG.Board.prototype.acos = Math.acos;
JXG.Board.prototype.asin = Math.asin;
JXG.Board.prototype.atan = Math.atan;
JXG.Board.prototype.ceil = Math.ceil;
JXG.Board.prototype.cos = Math.cos;
JXG.Board.prototype.exp = Math.exp;
JXG.Board.prototype.floor = Math.floor;
JXG.Board.prototype.log = Math.log;
JXG.Board.prototype.max = Math.max;
JXG.Board.prototype.min = Math.min;
JXG.Board.prototype.random = Math.random;
JXG.Board.prototype.sin = Math.sin;
JXG.Board.prototype.sqrt = Math.sqrt;
JXG.Board.prototype.tan = Math.tan;
JXG.Board.prototype.trunc = Math.ceil;

JXG.Board.prototype.factorial = function(n){ return JXG.Math.factorial(n); };
JXG.Board.prototype.binomial = function(n,k){ return JXG.Math.binomial(n,k); };

// Some shortcuts 
JXG.Point.prototype.setPositionX = function (method, x) {
    var y = (method==JXG.COORDS_BY_USER)?this.coords.usrCoords[2]:this.coords.scrCoords[2];
    this.setPosition(method,x,y);
};
JXG.Point.prototype.setPositionY = function (method, y) {
    var x = (method==JXG.COORDS_BY_USER)?this.coords.usrCoords[1]:this.coords.scrCoords[1];
    this.setPosition(method,x,y);
};
JXG.Board.prototype.getElement = function (el) {return JXG.getReference(this,el); };

/**
 * GUI interface
 **/
JXG.Board.prototype.intersectionOptions = ['point',[[JXG.OBJECT_CLASS_LINE,JXG.OBJECT_CLASS_LINE],[JXG.OBJECT_CLASS_LINE,JXG.OBJECT_CLASS_CIRCLE],[JXG.OBJECT_CLASS_CIRCLE,JXG.OBJECT_CLASS_CIRCLE]]];
JXG.Board.prototype.intersection = function(el1,el2,i,j){ 
    el1 = JXG.getReference(this,el1);
    el2 = JXG.getReference(this,el2);
    if (el1.elementClass==JXG.OBJECT_CLASS_CURVE || el2.elementClass==JXG.OBJECT_CLASS_CURVE) {
        return function(){return el1.board.algebra.meetCurveCurve(el1,el2,i,j); };
    } else {
        return function(){return el1.board.algebra.meet(el1.stdform,el2.stdform,i); };
    }
}; //returns a single point of intersection
JXG.Board.prototype.intersectionFunc = function(el1,el2,i,j){ return this.intersection(el1,el2,i,j); }; 

/**
* Intersectionof  circles and line
*/ 
JXG.Board.prototype.otherIntersection = function(el1,el2,el){ 
    el1 = JXG.getReference(this,el1);
    el2 = JXG.getReference(this,el2);
    return function(){
        var c = el1.board.algebra.meet(el1.stdform,el2.stdform,0);
        if (Math.abs(el.X()-c.usrCoords[1])>JXG.Math.eps ||
            Math.abs(el.Y()-c.usrCoords[2])>JXG.Math.eps ||
            Math.abs(el.Z()-c.usrCoords[0])>JXG.Math.eps) {
            return c;
        } else {
            return el1.board.algebra.meet(el1.stdform,el2.stdform,1);
        }
    };
}; //returns a single point of intersection


JXG.Board.prototype.pointFunc = function(){return [null];};
JXG.Board.prototype.pointOptions = ['point',[[JXG.OBJECT_CLASS_POINT]]];

JXG.Board.prototype.lineFunc = function(){return arguments;};
JXG.Board.prototype.lineOptions = ['line',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT]]];

JXG.Board.prototype.linesegmentFunc = function(){return arguments;};
JXG.Board.prototype.linesegmentOptions = ['line',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT]]];
JXG.Board.prototype.linesegmentAtts = {straightFirst : false, straightLast : false };

JXG.Board.prototype.arrowFunc = function(){return arguments;};
JXG.Board.prototype.arrowOptions = ['arrow',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT]]];

JXG.Board.prototype.circleFunc = function(){return arguments;};
JXG.Board.prototype.circleOptions = ['circle',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT],[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_LINE],[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_CIRCLE]]];

JXG.Board.prototype.arrowparallelOptions = ['arrowparallel',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_LINE]]];
JXG.Board.prototype.arrowparallelFunc = function(){return arguments;};

JXG.Board.prototype.bisectorOptions = ['bisector',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT]]];
JXG.Board.prototype.bisectorFunc = function(){return arguments;};

JXG.Board.prototype.circumcircleOptions = ['circumcircle',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT]]];
JXG.Board.prototype.circumcircleFunc = function(){return arguments;};

JXG.Board.prototype.circumcirclemidpointOptions = ['circumcirclemidpoint',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT]]];
JXG.Board.prototype.circumcirclemidpointFunc = function(){return arguments;};

JXG.Board.prototype.integralOptions = ['integral',[[]]];
JXG.Board.prototype.integralFunc = function(){return arguments;};

JXG.Board.prototype.midpointOptions = ['midpoint',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT],[JXG.OBJECT_CLASS_LINE]]];
JXG.Board.prototype.midpointFunc = function(){return arguments;};

JXG.Board.prototype.mirrorpointOptions = ['mirrorpoint',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT]]];
JXG.Board.prototype.mirrorpointFunc = function(){return arguments;};

JXG.Board.prototype.normalOptions = ['normal',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_LINE]]];
JXG.Board.prototype.normalFunc = function(){return arguments;};

JXG.Board.prototype.parallelOptions = ['parallel',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_LINE]]];
JXG.Board.prototype.parallelFunc = function(){return arguments;};

JXG.Board.prototype.parallelpointOptions = ['parallelpoint',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_POINT]]];
JXG.Board.prototype.parallelpointFunc = function(){return arguments;};

JXG.Board.prototype.perpendicularOptions = ['perpendicular',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_LINE]]];
JXG.Board.prototype.perpendicularFunc = function(){return arguments;};

JXG.Board.prototype.perpendicularpointOptions = ['perpendicularpoint',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_LINE]]];
JXG.Board.prototype.perpendicularpointFunc = function(){return arguments;};

JXG.Board.prototype.reflectionOptions = ['reflection',[[JXG.OBJECT_CLASS_POINT,JXG.OBJECT_CLASS_LINE]]];
JXG.Board.prototype.reflectionFunc = function(){return arguments;};

// Wrapper for not-singleton-pstricks. this could be removed after the next release
// and adjusting examples/pstricks.html and pstricks example in the wiki
// (http://jsxgraph.uni-bayreuth.de/wiki/index.php/PsTricks_export)
JXG.Board.prototype.pstricks = {};
JXG.Board.prototype.pstricks.givePsTricksToDiv = function(divId, board) {
    JXG.PsTricks.givePsTricksToDiv(divId, board);
};

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview In this file the geometry object Ticks is defined. Ticks provides
 * methods for creation and management of ticks on an axis.
 * @author graphjs
 * @version 0.1
 */

/**
 * Creates ticks for an axis.
 * @class Ticks provides methods for creation and management
 * of ticks on an axis.
 * @param {JXG.Line} line Reference to the axis the ticks are drawn on.
 * @param {Number,Array,Function} ticks Number, array or function defining the ticks.
 * @param {int} major Every major-th tick is drawn with heightmajorHeight, the other ones are drawn with height minorHeight.
 * @param {int} majorHeight The height used to draw major ticks.
 * @param {int} minorHeight The height used to draw minor ticks.
 * @param {String} id Unique identifier for this object.  If null or an empty string is given,
 * an unique id will be generated by Board.
 * @param {String} name Not necessarily unique name, won't be visible or used by this object.
 * @see JXG.Board#addTicks
 * @constructor
 * @extends JXG.GeometryElement
 */
JXG.Ticks = function (line, ticks, minor, majorHeight, minorHeight, id, name, layer) {
    /* Call the constructor of GeometryElement */
    this.constructor();

    /**
     * Type of GeometryElement, value is OBJECT_TYPE_ARC.
     * @final
     * @type int
     */
    this.type = JXG.OBJECT_TYPE_TICKS;

    /**
     * Class of the element, value is OBJECT_CLASS_CIRCLE.
     * @final
     * @type int
     */
    this.elementClass = JXG.OBJECT_CLASS_OTHER;

    /**
     * Set the display layer.
     */
    //if (layer == null) layer = board.options.layer['line']; // no board available
    //this.layer = layer;

    /**
     * The line the ticks belong to.
     * @type JXG.Line
     */
    this.line = line;

    /**
     * The board the ticks line is drawn on.
     * @type JXG.Board
     */
    this.board = this.line.board;

    /**
     * A function calculating ticks delta depending on the ticks number.
     * @type Function
     */
    this.ticksFunction = null;

    /**
     * Array of fixed ticks.
     * @type Array
     */
    this.fixedTicks = null;

    /**
     * Equidistant ticks. Distance is defined by ticksFunction
     * @type bool
     */
    this.equidistant = false;

    if(JXG.isFunction(ticks))
        this.ticksFunction = ticks;
    else if(JXG.isArray(ticks))
        this.fixedTicks = ticks;
    else {
        if(Math.abs(ticks) < JXG.Math.eps)
            ticks = this.board.options.line.ticks.defaultDistance;
        this.ticksFunction = function (i) { return ticks; };
        this.equidistant = true;
    }

    /**
     * minorTicks is the number of minor ticks between two major ticks.
     * @type int
     */
    this.minorTicks = ( (minor == null)? this.board.options.line.ticks.minorTicks : minor);
    if(this.minorTicks < 0)
        this.minorTicks = -this.minorTicks;

    /**
     * Total height of a major tick.
     * @type int
     */
    this.majorHeight = ( (majorHeight == null) || (majorHeight == 0) ? this.board.options.line.ticks.majorHeight : majorHeight);
    if(this.majorHeight < 0)
        this.majorHeight = -this.majorHeight;

    /**
     * Total height of a minor tick.
     * @type int
     */
    this.minorHeight = ( (minorHeight == null) || (minorHeight == 0) ? this.board.options.line.ticks.minorHeight : minorHeight);
    if(this.minorHeight < 0)
        this.minorHeight = -this.minorHeight;

    /**
     * Least distance between two ticks, measured in pixels.
     * @type int
     */
    this.minTicksDistance = this.board.options.line.ticks.minTicksDistance;

    /**
     * Maximum distance between two ticks, measured in pixels. Is used only when insertTicks
     * is set to true.
     * @type int
     * @see #insertTicks
     */
    this.maxTicksDistance = this.board.options.line.ticks.maxTicksDistance;

    /**
     * If the distance between two ticks is too big we could insert new ticks. If insertTicks
     * is <tt>true</tt>, we'll do so, otherwise we leave the distance as is.
     * This option is ignored if equidistant is false.
     * @type bool
     * @see #equidistant
     * @see #maxTicksDistance
     */
    this.insertTicks = this.board.options.line.ticks.insertTicks;

    /**
     * Draw the zero tick, that lies at line.point1?
     * @type bool
     */
    this.drawZero = this.board.options.line.ticks.drawZero;

    /**
     * Draw labels yes/no
     * @type bool
     */
    this.drawLabels = this.board.options.line.ticks.drawLabels;

    /**
     * Array where the labels are saved. There is an array element for every tick,
     * even for minor ticks which don't have labels. In this case the array element
     * contains just <tt>null</tt>.
     * @type array
     */
    this.labels = [];

    /* Call init defined in GeometryElement to set board, id and name property */
    this.init(this.board, id, name);

    this.visProp['visible'] = true;

    this.visProp['fillColor'] = this.line.visProp['fillColor'];
    this.visProp['highlightFillColor'] = this.line.visProp['highlightFillColor'];
    this.visProp['strokeColor'] = this.line.visProp['strokeColor'];
    this.visProp['highlightStrokeColor'] = this.line.visProp['highlightStrokeColor'];
    this.visProp['strokeWidth'] = this.line.visProp['strokeWidth'];

    /* Register ticks at line. */
    this.id = this.line.addTicks(this);
};

JXG.Ticks.prototype = new JXG.GeometryElement;

/**
 * Always returns false.
 * @param {int} x Coordinate in x direction, screen coordinates.
 * @param {int} y Coordinate in y direction, screen coordinates.
 * @return {bool} Always returns false.
 */
JXG.Ticks.prototype.hasPoint = function (x, y) {
   return false;
};

/**
 * This function acutally calculates the tick for one direction.
 * @param {JXG.Point} start Determines the start point from where the ticks are drawn.
 * This point is equal to the zero point of the line.
 * @param {JXG.Point} end Point to which the ticks are drawn.
 * @param {int} direction This determines the labels signum.
 * @param {int} over This should be 1 if we aren't allowed to draw ticks beyond the end point (if the line's a segment)
 * and should be zero if ticks should be drawn until the end of drawing area (if the line's a straight line).
 */
JXG.Ticks.prototype.makeTicks = function(start, end, direction, over) {
    // distance between start and end points
    var dx = start.usrCoords[1]-end.usrCoords[1]; // delta x
    var dy = start.usrCoords[2]-end.usrCoords[2]; // delta y

    // the distance between two ticks
    var ticksDelta = 0;

    // the length of the axis between start and end
    var total_length = Math.sqrt(dx*dx + dy*dy);

    if (total_length<=JXG.Math.eps)
        return;

    // x and y store the coordinates of the current tick to add
    var x = start.usrCoords[1];
    var y = start.usrCoords[2];

    // i is the amount of ticks drawn until now.
    var i = direction/Math.abs(direction);
    ticksDelta = Math.abs(this.ticksFunction(i));

    // this calculates the difference between the last and the current tick in
    // x and y direction in user coordinates
    var deltaX = (ticksDelta * dx) / (total_length);
    var deltaY = (ticksDelta * dy) / (total_length);
    
    // the current distance between two ticks in SCREEN coordinates
    // this is used to omit ticks, if they are too close to each other
    // and to add new ticks, if they are too far away from each other.
    var dist = 0;

    // if we have equidistant ticks do some delta-correction before we enter the main loop
    if(this.equidistant) {
        var deltaX_original = deltaX;
        var deltaY_original = deltaY;
        var ticksDelta_original = ticksDelta;
        
        var zero = new JXG.Coords(JXG.COORDS_BY_USER, [0, 0], this.board);
        var tmpCoords = new JXG.Coords(JXG.COORDS_BY_USER, [deltaX, deltaY], this.board);
        dist = (tmpCoords.scrCoords[1]-zero.scrCoords[1])*(tmpCoords.scrCoords[1]-zero.scrCoords[1]) + 
               (tmpCoords.scrCoords[2]-zero.scrCoords[2])*(tmpCoords.scrCoords[2]-zero.scrCoords[2]);

        ticksDelta = Math.pow(10,Math.floor(Math.log(ticksDelta)/Math.LN10));
        deltaX = (ticksDelta * dx) / (total_length);
        deltaY = (ticksDelta * dy) / (total_length);

        // If necessary, reduce ticksDelta
        while(dist > 8*this.minTicksDistance*this.minTicksDistance) {
            ticksDelta /= 10;
            deltaX = (ticksDelta * dx) / (total_length);
            deltaY = (ticksDelta * dy) / (total_length);

            tmpCoords = new JXG.Coords(JXG.COORDS_BY_USER, [deltaX, deltaY], this.board);
            dist = (tmpCoords.scrCoords[1]-zero.scrCoords[1])*(tmpCoords.scrCoords[1]-zero.scrCoords[1]) + 
                   (tmpCoords.scrCoords[2]-zero.scrCoords[2])*(tmpCoords.scrCoords[2]-zero.scrCoords[2]);
        }

        // If necessary, enlarge ticksDelta
        var factor = 5;
        while(dist < this.minTicksDistance*this.minTicksDistance) {
            ticksDelta *= factor;
            if (factor==5) { 
                factor = 2;
            } else {
                factor = 5;
            }
            deltaX = (ticksDelta * dx) / (total_length);
            deltaY = (ticksDelta * dy) / (total_length);

            tmpCoords = new JXG.Coords(JXG.COORDS_BY_USER, [deltaX, deltaY], this.board);
            dist = (tmpCoords.scrCoords[1]-zero.scrCoords[1])*(tmpCoords.scrCoords[1]-zero.scrCoords[1]) + 
                   (tmpCoords.scrCoords[2]-zero.scrCoords[2])*(tmpCoords.scrCoords[2]-zero.scrCoords[2]);
        }
    }

    // position is the current position on the axis
    var position = direction*ticksDelta;

    // reference to the last added tick coordinates object
    var lastTick = new JXG.Coords(JXG.COORDS_BY_USER, [x,y], this.board);

    var newTick = null;

    // minor ticks
    var minTick = null;

    var labelText = '';

    // temporary label object to create labels for ticks
    var label = null;

    // run the loop at least one time
    var first = true;

    // we also need access to the last ticks coordinates
    // store the old coordinates
    var lastX = x;
    var lastY = y;

    // loop abort criteria is:
    // our next tick is completely out of sight.
    while( first || 
           (
             this.board.sgn(deltaX)*(x-over*deltaX) >= this.board.sgn(deltaX)*end.usrCoords[1] && 
             this.board.sgn(deltaY)*(y-over*deltaY) >= this.board.sgn(deltaY)*end.usrCoords[2]
            ) 
         ) {
        // we're in it, so we have at least one tick and the deltaX/Y correction
        // for equidistant ticks.
        first = false;

        // calculate the new ticks coordinates
        x = x - deltaX;
        y = y - deltaY;

        // and put them into a coords object
        newTick = new JXG.Coords(JXG.COORDS_BY_USER, [x,y], this.board);

        // we need to calculate the distance. if we're in equidistant mode, we only need
        // to calculate it once, otherwise on every walk through the loop.
        if(!this.equidistant) {
            dist = (lastTick.scrCoords[1]-newTick.scrCoords[1])*(lastTick.scrCoords[1]-newTick.scrCoords[1]) +
                   (lastTick.scrCoords[2]-newTick.scrCoords[2])*(lastTick.scrCoords[2]-newTick.scrCoords[2]);
        }

        // if we're in equidistant mode and want to insert additional ticks automatically, whenever
        // the distance between two ticks is too big, we need to calculate the new deltaX and deltaY.
        if(this.insertTicks && this.equidistant && (dist > this.maxTicksDistance*this.maxTicksDistance)) {
            // dist is indeed the distance squared. repeat this, until we fall below the maxTicksDistance limit
            while (dist > this.maxTicksDistance*this.maxTicksDistance) {
                // half the distance
                deltaX *= 0.5;
                deltaY *= 0.5;
                ticksDelta *= 0.5;

                // move back towards the zeropoint
                x += deltaX;
                y += deltaY;

                position = position - direction*ticksDelta;

                // recalculate newTick coordinates and distance
                newTick = new JXG.Coords(JXG.COORDS_BY_USER, [x,y], this.board);
                dist = (lastTick.scrCoords[1]-newTick.scrCoords[1])*(lastTick.scrCoords[1]-newTick.scrCoords[1]) +
                       (lastTick.scrCoords[2]-newTick.scrCoords[2])*(lastTick.scrCoords[2]-newTick.scrCoords[2]);
            }
        }

        if(this.equidistant) {
            for(var z=1; z<this.minorTicks+1; z++) {
                minTick = new JXG.Coords(JXG.COORDS_BY_USER, [lastX - (deltaX*z)/(this.minorTicks+1) ,lastY - (deltaY*z)/(this.minorTicks+1)], this.board);
                minTick.major=false;
                this.ticks.push(minTick);
                this.labels.push(null);
            }
        }

        // if we're not below the minimum distance between two ticks, add the tick to our list
        if(this.equidistant || (dist > this.minTicksDistance*this.minTicksDistance)) {
            newTick.major = true;

            this.ticks.push(newTick);
            labelText = position.toString();
            if(labelText.length > 5)
                labelText = position.toPrecision(3).toString();
            label = new JXG.Text(this.board, labelText, null, [newTick.usrCoords[1], newTick.usrCoords[2]], this.id+i+"Label", null, null, true, this.board.options.text.defaultType);
            label.distanceX = 0;
            label.distanceY = -10;
            //label.setCoordinates(newTick);
            /*label.coords = new JXG.Coords(JXG.COORDS_BY_USER,
                                    [newTick.usrCoords[1]*1+label.distanceX/(this.board.stretchX),
                                     newTick.usrCoords[2]*1+label.distanceY/(this.board.stretchY)],
                                    this.board);*/
            label.setCoords(newTick.usrCoords[1]*1+label.distanceX/(this.board.stretchX), 
                            newTick.usrCoords[2]*1+label.distanceY/(this.board.stretchY));
            if (this.drawLabels) {
                label.visProp['visible'] = true; 
            }
            else {
                label.visProp['visible'] = false;
            }
            this.labels.push(label);
            // store the old coordinates
            lastX = x;
            lastY = y;

            lastTick = newTick;
        }

        i = i + direction*1;
        
        // if not equidistant, calculate the distance to the next tick
        if(!this.equidistant) {
            ticksDelta = Math.abs(this.ticksFunction(i));
        }
        
        // calculate the ticks label text data
        position = position + direction*ticksDelta;

        // recalculate new delta* only if needed
        if(!this.equidistant) {
            deltaX = (ticksDelta * dx) / (total_length);
            deltaY = (ticksDelta * dy) / (total_length);
        }
    }
};

/**
 * (Re-)calculates the ticks coordinates.
 */
JXG.Ticks.prototype.calculateTicksCoordinates = function() {
    // in which direction do we have to go?
    // DIR_PLUS is from line.p1 to line.p2, DIR_MINUS
    // is the other direction.
    var DIR_MINUS = 1;
    var DIR_PLUS = 2;
    var required = DIR_MINUS + DIR_PLUS;

    // calculate start (c1) and end (c2) points
    // copy existing lines point coordinates
    var c1 = new JXG.Coords(JXG.COORDS_BY_USER, [this.line.point1.coords.usrCoords[1], this.line.point1.coords.usrCoords[2]], this.board);
    var c2 = new JXG.Coords(JXG.COORDS_BY_USER, [this.line.point2.coords.usrCoords[1], this.line.point2.coords.usrCoords[2]], this.board);

    // now let the renderer calculate start and end point of the line on the board, i.e.
    // intersection points of line with the boards edges in the case that the line is a straight.
    this.board.renderer.calcStraight(this.line, c1, c2);
    
    // point1 is our reference point (where zero lies on the axis)
    var p1 = this.line.point1.coords;

    // first we have to look at some special cases, i.e. if our zero point on line
    // is outside the viewable area, there are some ticks we don't need to look at.
    if(this.board.renderer.isSameDirection(p1, c1, c2)) {
        if(this.board.renderer.isSameDirection(p1, this.line.point2.coords, c1)) {
            required = DIR_PLUS;
            if(p1.distance(JXG.COORDS_BY_USER, c1) > p1.distance(JXG.COORDS_BY_USER, c2))
                c2 = c1;
        } else {
            required = DIR_MINUS;
            if(p1.distance(JXG.COORDS_BY_USER, c1) < p1.distance(JXG.COORDS_BY_USER, c2))
                c1 = c2;
        }
    } else /* p1 can be seen on the drawing area */ {
        // make sure, point2 is on that part of the line with positive tick markers on it
        if(this.board.renderer.isSameDirection(p1, this.line.point2.coords, c1)) {
            var ct = c1;
            c1 = c2;
            c2 = ct;
        }
    }
    
    if(this.ticks != null) {
        for(var j=0; j<this.ticks.length; j++) {
            if(this.labels[j] != null) {
                if (this.labels[j].visProp['visible']) this.board.renderer.remove(this.labels[j].rendNode);
            }
        }
    }

    // initialize storage arrays
    // ticks stores the ticks coordinates
    this.ticks = new Array();
    // labels stores the text to display beside the ticks
    this.labels = new Array();

    var label = null;
    var labelText = '';

    // reference to the newly added tick coordinates object
    var newTick = null;
    if(this.ticksFunction != null) {
        // add tick at p1?
        if(this.drawZero) {
            newTick = new JXG.Coords(JXG.COORDS_BY_USER, [p1.usrCoords[1], p1.usrCoords[2]], this.board);
            this.ticks.push(newTick);
            //label = new JXG.Label(this.board, "0", newTick, this.id+"0Label");
            label = new JXG.Text(this.board, "0", null, [p1.usrCoords[1], p1.usrCoords[2]], this.id+"0Label", null, null, true, this.board.options.text.defaultType);
            if (this.drawLabels) {
                label.visProp['visible'] = true; 
            }
            else {
                label.visProp['visible'] = false;
            }
            this.labels.push(label);

            this.ticks[0].major = true;
        }

        if(DIR_MINUS == (required & DIR_MINUS)) {
            if(this.line.visProp['straightFirst'])
                this.makeTicks(p1, c1, -1, 0);
// the negative part even doesn't exist on this line...
//            else
//                this.makeTicks(p1, c1, -1, 1);
        }
        if(DIR_PLUS == (required & DIR_PLUS)) {
            if(this.line.visProp['straightLast'])
                this.makeTicks(p1, c2, +1, 0);
            else {
                this.makeTicks(p1, this.line.point2.coords, +1, 1);
            }
        }
    } else {
        // we have an array of fixed ticks we have to draw
        if(!this.line.visProp['straightFirst'])
            c1 = p1;
        var dx_minus = p1.usrCoords[1]-c1.usrCoords[1];
        var dy_minus = p1.usrCoords[2]-c1.usrCoords[2];
        var length_minus = Math.sqrt(dx_minus*dx_minus + dy_minus*dy_minus);

        if(!this.line.visProp['straightLast'])
            c2 = this.line.point2.coords;
        var dx_plus = p1.usrCoords[1]-c2.usrCoords[1];
        var dy_plus = p1.usrCoords[2]-c2.usrCoords[2];
        var length_plus = Math.sqrt(dx_plus*dx_plus + dy_plus*dy_plus);

        // new ticks coordinates
        var nx = 0;
        var ny = 0;

        for(var i=0; i<this.fixedTicks.length; i++) {
            // is this tick visible?
            if((-length_minus <= this.fixedTicks[i]) && (this.fixedTicks[i] <= length_plus)) {
                if(this.fixedTicks[i] < 0) {
                    nx = Math.abs(dx_minus) * this.fixedTicks[i]/length_minus;
                    ny = Math.abs(dy_minus) * this.fixedTicks[i]/length_minus;
                } else {
                    nx = Math.abs(dx_plus) * this.fixedTicks[i]/length_plus;
                    ny = Math.abs(dy_plus) * this.fixedTicks[i]/length_plus;
                }

                newTick = new JXG.Coords(JXG.COORDS_BY_USER, [p1.usrCoords[1] + nx, p1.usrCoords[2] + ny], this.board);
                this.ticks.push(newTick);
                this.ticks[this.ticks.length-1].major = true;
                labelText = this.fixedTicks[i].toString();
                if(labelText.length > 5)
                    labelText = this.fixedTicks[i].toFixed(3).toString();
                //label = new JXG.Label(this.board, labelText, newTick, this.id+i+"Label");
                label = new JXG.Text(this.board, labelText, null, [p1.usrCoords[1] + nx, p1.usrCoords[2] + ny], this.id+i+"Label", null, null, true, this.board.options.text.defaultType);
                label.distanceX = 0;
                label.distanceY = -10;
                //label.setCoordinates(newTick);
                /*label.coords = new JXG.Coords(JXG.COORDS_BY_USER,
                             [newTick.usrCoords[1]*1+label.distanceX/(this.board.stretchX),
                              newTick.usrCoords[2]*1+label.distanceY/(this.board.stretchY)],
                             this.board);*/    
                label.setCoords(newTick.usrCoords[1]*1+label.distanceX/(this.board.stretchX), 
                                newTick.usrCoords[2]*1+label.distanceY/(this.board.stretchY));
                if (this.drawLabels) {
                    label.visProp['visible'] = true; 
                } else {
                    label.visProp['visible'] = false;
                }
                this.labels.push(label);
            }
        }
    }

    // this piece of code was in AbstractRenderer.updateAxisTicksInnerLoop
    // and has been moved in here to clean up the code.
    //
    // The code above only calculates the position of the ticks. The following code parts
    // calculate the dx and dy values which make ticks out of this positions, i.e. from the
    // position (p_x, p_y) calculated above we have to draw a line from
    // (p_x - dx, py - dy) to (p_x + dx, p_y + dy) to get a tick.
    var eps = JXG.Math.eps;
    var slope = -this.line.getSlope();

    var distMaj = this.majorHeight/2;
    var distMin = this.minorHeight/2;
    var dxMaj = 0; var dyMaj = 0;
    var dxMin = 0; var dyMin = 0;

    if(Math.abs(slope) < eps) {
        // if the slope of the line is (almost) 0, we can set dx and dy directly
        dxMaj = 0;
        dyMaj = distMaj;
        dxMin = 0;
        dyMin = distMin;
    } else if((Math.abs(slope) > 1/eps) || (isNaN(slope))) {
        // if the slope of the line is (theoretically) infinite, we can set dx and dy directly
        dxMaj = distMaj;
        dyMaj = 0;
        dxMin = distMin;
        dyMin = 0;
    } else {
        // here we have to calculate dx and dy depending on the slope and the length of the tick (dist)
        // if slope is the line's slope, the tick's slope is given by
        //
        //            1          dy
        //     -   -------  =   ----                 (I)
        //          slope        dx
        //
        // when dist is the length of the tick, using the pythagorean theorem we get
        //
        //     dx*dx + dy*dy = dist*dist             (II)
        //
        // dissolving (I) by dy and applying that to equation (II) we get the following formulas for dx and dy
        dxMaj = -distMaj/Math.sqrt(1/(slope*slope) + 1);
        dyMaj = dxMaj/slope;
        dxMin = -distMin/Math.sqrt(1/(slope*slope) + 1);
        dyMin = dxMin/slope;
    }
    this.board.renderer.updateTicks(this,dxMaj,dyMaj,dxMin,dyMin);
};

/**
 * Uses the boards renderer to update the arc.
 * update() is not needed for arc.
 */
JXG.Ticks.prototype.updateRenderer = function () {
    if (this.needsUpdate) {
        this.calculateTicksCoordinates();
        this.needsUpdate = false;
    }
};

/**
 * Creates new ticks.
 * @param {JXG.Board} board The board the ticks are put on.
 * @param {Array} parents Array containing a line and an array of positions, where ticks should be put on that line or
 *   a function that calculates the distance based on the ticks number that is given as a parameter. E.g.:<br />
 *   <tt>var ticksFunc = function(i) {</tt><br />
 *   <tt>    return 2;</tt><br />
 *   <tt>}</tt><br />
 *   for ticks with distance 2 between each tick.
 * @param {Object} attributs Object containing properties for the element such as stroke-color and visibility. See @see JXG.GeometryElement#setProperty
 * @type JXG.Ticks
 * @return Reference to the created ticks object.
 */
JXG.createTicks = function(board, parents, attributes) {
    var el;
    attributes = JXG.checkAttributes(attributes,{layer:null});
    if ( (parents[0].elementClass == JXG.OBJECT_CLASS_LINE) && (JXG.isFunction(parents[1]) || JXG.isArray(parents[1]) || JXG.isNumber(parents[1]))) {
        el = new JXG.Ticks(parents[0], parents[1], attributes['minorTicks'], attributes['majHeight'], attributes['minHeight'], attributes['id'], attributes['name'], attributes['layer']);
    } else
        throw new Error("JSXGraph: Can't create Ticks with parent types '" + (typeof parents[0]) + "' and '" + (typeof parents[1]) + "' and '" + (typeof parents[2]) + "'.");

    return el;
};

JXG.JSXGraph.registerElement('ticks', JXG.createTicks);

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview Utilities for uncompressing and base64 decoding
 */

/**
  * @class Util class
  * Class for gunzipping, unzipping and base64 decoding of files.
  * It is used for reading GEONExT, Geogebra and Intergeo files.
  *
  * Only Huffman codes are decoded in gunzip.
  * The code is based on the source code for gunzip.c by Pasi Ojala 
  * @see <a href="http://www.cs.tut.fi/~albert/Dev/gunzip/gunzip.c">http://www.cs.tut.fi/~albert/Dev/gunzip/gunzip.c</a>
  * @see <a href="http://www.cs.tut.fi/~albert">http://www.cs.tut.fi/~albert</a>
  */
JXG.Util = {};
                                 
/**
 * Unzip zip files
 */
JXG.Util.Unzip = function (barray){
    var outputArr = [],
        output = "",
        debug = false,
        gpflags,
        files = 0,
        unzipped = [],
        crc,
        buf32k = new Array(32768),
        bIdx = 0,
        modeZIP=false,

        CRC, SIZE,
    
        bitReverse = [
        0x00, 0x80, 0x40, 0xc0, 0x20, 0xa0, 0x60, 0xe0,
        0x10, 0x90, 0x50, 0xd0, 0x30, 0xb0, 0x70, 0xf0,
        0x08, 0x88, 0x48, 0xc8, 0x28, 0xa8, 0x68, 0xe8,
        0x18, 0x98, 0x58, 0xd8, 0x38, 0xb8, 0x78, 0xf8,
        0x04, 0x84, 0x44, 0xc4, 0x24, 0xa4, 0x64, 0xe4,
        0x14, 0x94, 0x54, 0xd4, 0x34, 0xb4, 0x74, 0xf4,
        0x0c, 0x8c, 0x4c, 0xcc, 0x2c, 0xac, 0x6c, 0xec,
        0x1c, 0x9c, 0x5c, 0xdc, 0x3c, 0xbc, 0x7c, 0xfc,
        0x02, 0x82, 0x42, 0xc2, 0x22, 0xa2, 0x62, 0xe2,
        0x12, 0x92, 0x52, 0xd2, 0x32, 0xb2, 0x72, 0xf2,
        0x0a, 0x8a, 0x4a, 0xca, 0x2a, 0xaa, 0x6a, 0xea,
        0x1a, 0x9a, 0x5a, 0xda, 0x3a, 0xba, 0x7a, 0xfa,
        0x06, 0x86, 0x46, 0xc6, 0x26, 0xa6, 0x66, 0xe6,
        0x16, 0x96, 0x56, 0xd6, 0x36, 0xb6, 0x76, 0xf6,
        0x0e, 0x8e, 0x4e, 0xce, 0x2e, 0xae, 0x6e, 0xee,
        0x1e, 0x9e, 0x5e, 0xde, 0x3e, 0xbe, 0x7e, 0xfe,
        0x01, 0x81, 0x41, 0xc1, 0x21, 0xa1, 0x61, 0xe1,
        0x11, 0x91, 0x51, 0xd1, 0x31, 0xb1, 0x71, 0xf1,
        0x09, 0x89, 0x49, 0xc9, 0x29, 0xa9, 0x69, 0xe9,
        0x19, 0x99, 0x59, 0xd9, 0x39, 0xb9, 0x79, 0xf9,
        0x05, 0x85, 0x45, 0xc5, 0x25, 0xa5, 0x65, 0xe5,
        0x15, 0x95, 0x55, 0xd5, 0x35, 0xb5, 0x75, 0xf5,
        0x0d, 0x8d, 0x4d, 0xcd, 0x2d, 0xad, 0x6d, 0xed,
        0x1d, 0x9d, 0x5d, 0xdd, 0x3d, 0xbd, 0x7d, 0xfd,
        0x03, 0x83, 0x43, 0xc3, 0x23, 0xa3, 0x63, 0xe3,
        0x13, 0x93, 0x53, 0xd3, 0x33, 0xb3, 0x73, 0xf3,
        0x0b, 0x8b, 0x4b, 0xcb, 0x2b, 0xab, 0x6b, 0xeb,
        0x1b, 0x9b, 0x5b, 0xdb, 0x3b, 0xbb, 0x7b, 0xfb,
        0x07, 0x87, 0x47, 0xc7, 0x27, 0xa7, 0x67, 0xe7,
        0x17, 0x97, 0x57, 0xd7, 0x37, 0xb7, 0x77, 0xf7,
        0x0f, 0x8f, 0x4f, 0xcf, 0x2f, 0xaf, 0x6f, 0xef,
        0x1f, 0x9f, 0x5f, 0xdf, 0x3f, 0xbf, 0x7f, 0xff
    ],
    
    cplens = [
        3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31,
        35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258, 0, 0
    ],

    cplext = [
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0, 99, 99
    ], /* 99==invalid */

    cpdist = [
        0x0001, 0x0002, 0x0003, 0x0004, 0x0005, 0x0007, 0x0009, 0x000d,
        0x0011, 0x0019, 0x0021, 0x0031, 0x0041, 0x0061, 0x0081, 0x00c1,
        0x0101, 0x0181, 0x0201, 0x0301, 0x0401, 0x0601, 0x0801, 0x0c01,
        0x1001, 0x1801, 0x2001, 0x3001, 0x4001, 0x6001
    ],

    cpdext = [
        0,  0,  0,  0,  1,  1,  2,  2,
        3,  3,  4,  4,  5,  5,  6,  6,
        7,  7,  8,  8,  9,  9, 10, 10,
        11, 11, 12, 12, 13, 13
    ],
    
    border = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15],
    
    bA = barray,

    bytepos=0,
    bitpos=0,
    bb = 1,
    bits=0,
    
    NAMEMAX = 256,
    
    nameBuf = [],
    
    fileout;
    
    function readByte(){
        var len = bA.length;
        bits+=8;
        if (bytepos<len){
            if (debug)
                document.write(bytepos+": "+bA[bytepos]+"<br>");
            return bA[bytepos++];
        } else
            return -1;
    };

    function byteAlign(){
        bb = 1;
    };
    
    function readBit(){
        var carry;
        bits++;
        carry = (bb & 1);
        bb >>= 1;
        if (bb==0){
            bb = readByte();
            carry = (bb & 1);
            bb = (bb>>1) | 0x80;
        }
        return carry;
    };

    function readBits(a) {
        var res = 0,
            i = a;
    
        while(i--) {
            res = (res<<1) | readBit();
        }
        if(a) {
            res = bitReverse[res]>>(8-a);
        }
        return res;
    };
        
    function flushBuffer(){
        //document.write('FLUSHBUFFER:'+buf32k);
        bIdx = 0;
    };
    function addBuffer(a){
        SIZE++;
        //CRC=updcrc(a,crc);
        buf32k[bIdx++] = a;
        outputArr.push(String.fromCharCode(a));
        //output+=String.fromCharCode(a);
        if(bIdx==0x8000){
            //document.write('ADDBUFFER:'+buf32k);
            bIdx=0;
        }
    };
    
    function HufNode() {
        this.b0=0;
        this.b1=0;
        this.jump = null;
        this.jumppos = -1;
    };

    var LITERALS = 288;
    
    var literalTree = new Array(LITERALS);
    var distanceTree = new Array(32);
    var treepos=0;
    var Places = null;
    var Places2 = null;
    
    var impDistanceTree = new Array(64);
    var impLengthTree = new Array(64);
    
    var len = 0;
    var fpos = new Array(17);
    fpos[0]=0;
    var flens;
    var fmax;
    
    function IsPat() {
        while (1) {
            if (fpos[len] >= fmax)
                return -1;
            if (flens[fpos[len]] == len)
                return fpos[len]++;
            fpos[len]++;
        }
    };

    function Rec() {
        var curplace = Places[treepos];
        var tmp;
        if (debug)
    		document.write("<br>len:"+len+" treepos:"+treepos);
        if(len==17) { //war 17
            return -1;
        }
        treepos++;
        len++;
    	
        tmp = IsPat();
        if (debug)
        	document.write("<br>IsPat "+tmp);
        if(tmp >= 0) {
            curplace.b0 = tmp;    /* leaf cell for 0-bit */
            if (debug)
            	document.write("<br>b0 "+curplace.b0);
        } else {
        /* Not a Leaf cell */
        curplace.b0 = 0x8000;
        if (debug)
        	document.write("<br>b0 "+curplace.b0);
        if(Rec())
            return -1;
        }
        tmp = IsPat();
        if(tmp >= 0) {
            curplace.b1 = tmp;    /* leaf cell for 1-bit */
            if (debug)
            	document.write("<br>b1 "+curplace.b1);
            curplace.jump = null;    /* Just for the display routine */
        } else {
            /* Not a Leaf cell */
            curplace.b1 = 0x8000;
            if (debug)
            	document.write("<br>b1 "+curplace.b1);
            curplace.jump = Places[treepos];
            curplace.jumppos = treepos;
            if(Rec())
                return -1;
        }
        len--;
        return 0;
    };

    function CreateTree(currentTree, numval, lengths, show) {
        var i;
        /* Create the Huffman decode tree/table */
        //document.write("<br>createtree<br>");
        if (debug)
        	document.write("currentTree "+currentTree+" numval "+numval+" lengths "+lengths+" show "+show);
        Places = currentTree;
        treepos=0;
        flens = lengths;
        fmax  = numval;
        for (i=0;i<17;i++)
            fpos[i] = 0;
        len = 0;
        if(Rec()) {
            //fprintf(stderr, "invalid huffman tree\n");
            if (debug)
            	alert("invalid huffman tree\n");
            return -1;
        }
        if (debug){
        	document.write('<br>Tree: '+Places.length);
        	for (var a=0;a<32;a++){
            	document.write("Places["+a+"].b0="+Places[a].b0+"<br>");
            	document.write("Places["+a+"].b1="+Places[a].b1+"<br>");
        	}
        }
    
        /*if(show) {
            var tmp;
            for(tmp=currentTree;tmp<Places;tmp++) {
                fprintf(stdout, "0x%03x  0x%03x (0x%04x)",tmp-currentTree, tmp->jump?tmp->jump-currentTree:0,(tmp->jump?tmp->jump-currentTree:0)*6+0xcf0);
                if(!(tmp.b0 & 0x8000)) {
                    //fprintf(stdout, "  0x%03x (%c)", tmp->b0,(tmp->b0<256 && isprint(tmp->b0))?tmp->b0:'�');
                }
                if(!(tmp.b1 & 0x8000)) {
                    if((tmp.b0 & 0x8000))
                        fprintf(stdout, "           ");
                    fprintf(stdout, "  0x%03x (%c)", tmp->b1,(tmp->b1<256 && isprint(tmp->b1))?tmp->b1:'�');
                }
                fprintf(stdout, "\n");
            }
        }*/
        return 0;
    };
    
    function DecodeValue(currentTree) {
        var len, i,
            xtreepos=0,
            X = currentTree[xtreepos],
            b;

        /* decode one symbol of the data */
        while(1) {
            b=readBit();
            if (debug)
            	document.write("b="+b);
            if(b) {
                if(!(X.b1 & 0x8000)){
                	if (debug)
                    	document.write("ret1");
                    return X.b1;    /* If leaf node, return data */
                }
                X = X.jump;
                len = currentTree.length;
                for (i=0;i<len;i++){
                    if (currentTree[i]===X){
                        xtreepos=i;
                        break;
                    }
                }
                //xtreepos++;
            } else {
                if(!(X.b0 & 0x8000)){
                	if (debug)
                    	document.write("ret2");
                    return X.b0;    /* If leaf node, return data */
                }
                //X++; //??????????????????
                xtreepos++;
                X = currentTree[xtreepos];
            }
        }
        if (debug)
        	document.write("ret3");
        return -1;
    };
    
    function DeflateLoop() {
    var last, c, type, i, len;

    do {
        /*if((last = readBit())){
            fprintf(errfp, "Last Block: ");
        } else {
            fprintf(errfp, "Not Last Block: ");
        }*/
        last = readBit();
        type = readBits(2);
        switch(type) {
            case 0:
            	if (debug)
                	alert("Stored\n");
                break;
            case 1:
            	if (debug)
                	alert("Fixed Huffman codes\n");
                break;
            case 2:
            	if (debug)
                	alert("Dynamic Huffman codes\n");
                break;
            case 3:
            	if (debug)
                	alert("Reserved block type!!\n");
                break;
            default:
            	if (debug)
                	alert("Unexpected value %d!\n", type);
                break;
        }

        if(type==0) {
            var blockLen, cSum;

            // Stored 
            byteAlign();
            blockLen = readByte();
            blockLen |= (readByte()<<8);

            cSum = readByte();
            cSum |= (readByte()<<8);

            if(((blockLen ^ ~cSum) & 0xffff)) {
                document.write("BlockLen checksum mismatch\n");
            }
            while(blockLen--) {
                c = readByte();
                addBuffer(c);
            }
        } else if(type==1) {
            var j;

            /* Fixed Huffman tables -- fixed decode routine */
            while(1) {
            /*
                256    0000000        0
                :   :     :
                279    0010111        23
                0   00110000    48
                :    :      :
                143    10111111    191
                280 11000000    192
                :    :      :
                287 11000111    199
                144    110010000    400
                :    :       :
                255    111111111    511
    
                Note the bit order!
                */

            j = (bitReverse[readBits(7)]>>1);
            if(j > 23) {
                j = (j<<1) | readBit();    /* 48..255 */

                if(j > 199) {    /* 200..255 */
                    j -= 128;    /*  72..127 */
                    j = (j<<1) | readBit();        /* 144..255 << */
                } else {        /*  48..199 */
                    j -= 48;    /*   0..151 */
                    if(j > 143) {
                        j = j+136;    /* 280..287 << */
                        /*   0..143 << */
                    }
                }
            } else {    /*   0..23 */
                j += 256;    /* 256..279 << */
            }
            if(j < 256) {
                addBuffer(j);
                //document.write("out:"+String.fromCharCode(j));
                /*fprintf(errfp, "@%d %02x\n", SIZE, j);*/
            } else if(j == 256) {
                /* EOF */
                break;
            } else {
                var len, dist;

                j -= 256 + 1;    /* bytes + EOF */
                len = readBits(cplext[j]) + cplens[j];

                j = bitReverse[readBits(5)]>>3;
                if(cpdext[j] > 8) {
                    dist = readBits(8);
                    dist |= (readBits(cpdext[j]-8)<<8);
                } else {
                    dist = readBits(cpdext[j]);
                }
                dist += cpdist[j];

                /*fprintf(errfp, "@%d (l%02x,d%04x)\n", SIZE, len, dist);*/
                for(j=0;j<len;j++) {
                    var c = buf32k[(bIdx - dist) & 0x7fff];
                    addBuffer(c);
                }
            }
            } // while
        } else if(type==2) {
            var j, n, literalCodes, distCodes, lenCodes;
            var ll = new Array(288+32);    // "static" just to preserve stack
    
            // Dynamic Huffman tables 
    
            literalCodes = 257 + readBits(5);
            distCodes = 1 + readBits(5);
            lenCodes = 4 + readBits(4);
            //document.write("<br>param: "+literalCodes+" "+distCodes+" "+lenCodes+"<br>");
            for(j=0; j<19; j++) {
                ll[j] = 0;
            }
    
            // Get the decode tree code lengths
    
            //document.write("<br>");
            for(j=0; j<lenCodes; j++) {
                ll[border[j]] = readBits(3);
                //document.write(ll[border[j]]+" ");
            }
            //fprintf(errfp, "\n");
            //document.write('<br>ll:'+ll);
            len = distanceTree.length;
            for (i=0; i<len; i++)
                distanceTree[i]=new HufNode();
            if(CreateTree(distanceTree, 19, ll, 0)) {
                flushBuffer();
                return 1;
            }
            if (debug){
            	document.write("<br>distanceTree");
            	for(var a=0;a<distanceTree.length;a++){
                	document.write("<br>"+distanceTree[a].b0+" "+distanceTree[a].b1+" "+distanceTree[a].jump+" "+distanceTree[a].jumppos);
                	/*if (distanceTree[a].jumppos!=-1)
                    	document.write(" "+distanceTree[a].jump.b0+" "+distanceTree[a].jump.b1);
                	*/
            	}
            }
            //document.write('<BR>tree created');
    
            //read in literal and distance code lengths
            n = literalCodes + distCodes;
            i = 0;
            var z=-1;
            if (debug)
            	document.write("<br>n="+n+" bits: "+bits+"<br>");
            while(i < n) {
                z++;
                j = DecodeValue(distanceTree);
                if (debug)
                	document.write("<br>"+z+" i:"+i+" decode: "+j+"    bits "+bits+"<br>");
                if(j<16) {    // length of code in bits (0..15)
                       ll[i++] = j;
                } else if(j==16) {    // repeat last length 3 to 6 times 
                       var l;
                    j = 3 + readBits(2);
                    if(i+j > n) {
                        flushBuffer();
                        return 1;
                    }
                    l = i ? ll[i-1] : 0;
                    while(j--) {
                        ll[i++] = l;
                    }
                } else {
                    if(j==17) {        // 3 to 10 zero length codes
                        j = 3 + readBits(3);
                    } else {        // j == 18: 11 to 138 zero length codes 
                        j = 11 + readBits(7);
                    }
                    if(i+j > n) {
                        flushBuffer();
                        return 1;
                    }
                    while(j--) {
                        ll[i++] = 0;
                    }
                }
            }
            /*for(j=0; j<literalCodes+distCodes; j++) {
                //fprintf(errfp, "%d ", ll[j]);
                if ((j&7)==7)
                    fprintf(errfp, "\n");
            }
            fprintf(errfp, "\n");*/
            // Can overwrite tree decode tree as it is not used anymore
            len = literalTree.length;
            for (i=0; i<len; i++)
                literalTree[i]=new HufNode();
            if(CreateTree(literalTree, literalCodes, ll, 0)) {
                flushBuffer();
                return 1;
            }
            len = literalTree.length;
            for (i=0; i<len; i++)
                distanceTree[i]=new HufNode();
            var ll2 = new Array();
            for (i=literalCodes; i <ll.length; i++){
                ll2[i-literalCodes]=ll[i];
            }    
            if(CreateTree(distanceTree, distCodes, ll2, 0)) {
                flushBuffer();
                return 1;
            }
            if (debug)
           		document.write("<br>literalTree");
            while(1) {
                j = DecodeValue(literalTree);
                if(j >= 256) {        // In C64: if carry set
                    var len, dist;
                    j -= 256;
                    if(j == 0) {
                        // EOF
                        break;
                    }
                    j--;
                    len = readBits(cplext[j]) + cplens[j];
    
                    j = DecodeValue(distanceTree);
                    if(cpdext[j] > 8) {
                        dist = readBits(8);
                        dist |= (readBits(cpdext[j]-8)<<8);
                    } else {
                        dist = readBits(cpdext[j]);
                    }
                    dist += cpdist[j];
                    while(len--) {
                        var c = buf32k[(bIdx - dist) & 0x7fff];
                        addBuffer(c);
                    }
                } else {
                    addBuffer(j);
                }
            }
        }
    } while(!last);
    flushBuffer();

    byteAlign();
    return 0;
};

JXG.Util.Unzip.prototype.unzipFile = function(name) {
    var i;
	this.unzip();
	//alert(unzipped[0][1]);
	for (i=0;i<unzipped.length;i++){
		if(unzipped[i][1]==name) {
			return unzipped[i][0];
		}
	}
	
  };
    
    
JXG.Util.Unzip.prototype.unzip = function() {
	//convertToByteArray(input);
	if (debug)
		alert(bA);
	/*for (i=0;i<bA.length*8;i++){
		document.write(readBit());
		if ((i+1)%8==0)
			document.write(" ");
	}*/
	/*for (i=0;i<bA.length;i++){
		document.write(readByte()+" ");
		if ((i+1)%8==0)
			document.write(" ");
	}
	for (i=0;i<bA.length;i++){
		document.write(bA[i]+" ");
		if ((i+1)%16==0)
			document.write("<br>");
	}	
	*/
	//alert(bA);
	nextFile();
	return unzipped;
  };
    
 function nextFile(){
 	if (debug)
 		alert("NEXTFILE");
 	outputArr = [];
 	var tmp = [];
 	modeZIP = false;
	tmp[0] = readByte();
	tmp[1] = readByte();
	if (debug)
		alert("type: "+tmp[0]+" "+tmp[1]);
	if (tmp[0] == parseInt("78",16) && tmp[1] == parseInt("da",16)){ //GZIP
		if (debug)
			alert("GEONExT-GZIP");
		DeflateLoop();
		if (debug)
			alert(outputArr.join(''));
		unzipped[files] = new Array(2);
    	unzipped[files][0] = outputArr.join('');
    	unzipped[files][1] = "geonext.gxt";
    	files++;
	}
	if (tmp[0] == parseInt("50",16) && tmp[1] == parseInt("4b",16)){ //ZIP
		modeZIP = true;
		tmp[2] = readByte();
		tmp[3] = readByte();
		if (tmp[2] == parseInt("3",16) && tmp[3] == parseInt("4",16)){
			//MODE_ZIP
			tmp[0] = readByte();
			tmp[1] = readByte();
			if (debug)
				alert("ZIP-Version: "+tmp[1]+" "+tmp[0]/10+"."+tmp[0]%10);
			
			gpflags = readByte();
			gpflags |= (readByte()<<8);
			if (debug)
				alert("gpflags: "+gpflags);
			
			var method = readByte();
			method |= (readByte()<<8);
			if (debug)
				alert("method: "+method);
			
			readByte();
			readByte();
			readByte();
			readByte();
			
			var crc = readByte();
			crc |= (readByte()<<8);
			crc |= (readByte()<<16);
			crc |= (readByte()<<24);
			
			var compSize = readByte();
			compSize |= (readByte()<<8);
			compSize |= (readByte()<<16);
			compSize |= (readByte()<<24);
			
			var size = readByte();
			size |= (readByte()<<8);
			size |= (readByte()<<16);
			size |= (readByte()<<24);
			
			if (debug)
				alert("local CRC: "+crc+"\nlocal Size: "+size+"\nlocal CompSize: "+compSize);
			
			var filelen = readByte();
			filelen |= (readByte()<<8);
			
			var extralen = readByte();
			extralen |= (readByte()<<8);
			
			if (debug)
				alert("filelen "+filelen);
			i = 0;
			nameBuf = [];
			while (filelen--){ 
				var c = readByte();
				if (c == "/" | c ==":"){
					i = 0;
				} else if (i < NAMEMAX-1)
					nameBuf[i++] = String.fromCharCode(c);
			}
			if (debug)
				alert("nameBuf: "+nameBuf);
			
			//nameBuf[i] = "\0";
			if (!fileout)
				fileout = nameBuf;
			
			var i = 0;
			while (i < extralen){
				c = readByte();
				i++;
			}
				
			CRC = 0xffffffff;
			SIZE = 0;
			
			if (size = 0 && fileOut.charAt(fileout.length-1)=="/"){
				//skipdir
				if (debug)
					alert("skipdir");
			}
			if (method == 8){
				DeflateLoop();
				if (debug)
					alert(outputArr.join(''));
				unzipped[files] = new Array(2);
				unzipped[files][0] = outputArr.join('');
    			unzipped[files][1] = nameBuf.join('');
    			files++;
				//return outputArr.join('');
			}
			skipdir();
		}
	}
 };
	
function skipdir(){
    var crc, 
        tmp = [],
        compSize, size, os, i, c;
    
	if ((gpflags & 8)) {
		tmp[0] = readByte();
		tmp[1] = readByte();
		tmp[2] = readByte();
		tmp[3] = readByte();
		
		if (tmp[0] == parseInt("50",16) && 
            tmp[1] == parseInt("4b",16) && 
            tmp[2] == parseInt("07",16) && 
            tmp[3] == parseInt("08",16))
        {
            crc = readByte();
            crc |= (readByte()<<8);
            crc |= (readByte()<<16);
            crc |= (readByte()<<24);
		} else {
			crc = tmp[0] | (tmp[1]<<8) | (tmp[2]<<16) | (tmp[3]<<24);
		}
		
		compSize = readByte();
		compSize |= (readByte()<<8);
		compSize |= (readByte()<<16);
		compSize |= (readByte()<<24);
		
		size = readByte();
		size |= (readByte()<<8);
		size |= (readByte()<<16);
		size |= (readByte()<<24);
		
		if (debug)
			alert("CRC:");
	}

	if (modeZIP)
		nextFile();
	
	tmp[0] = readByte();
	if (tmp[0] != 8) {
		if (debug)
			alert("Unknown compression method!");
        return 0;	
	}
	
	gpflags = readByte();
	if (debug){
		if ((gpflags & ~(parseInt("1f",16))))
			alert("Unknown flags set!");
	}
	
	readByte();
	readByte();
	readByte();
	readByte();
	
	readByte();
	os = readByte();
	
	if ((gpflags & 4)){
		tmp[0] = readByte();
		tmp[2] = readByte();
		len = tmp[0] + 256*tmp[1];
		if (debug)
			alert("Extra field size: "+len);
		for (i=0;i<len;i++)
			readByte();
	}
	
	if ((gpflags & 8)){
		i=0;
		nameBuf=[];
		while (c=readByte()){
			if(c == "7" || c == ":")
				i=0;
			if (i<NAMEMAX-1)
				nameBuf[i++] = c;
		}
		//nameBuf[i] = "\0";
		if (debug)
			alert("original file name: "+nameBuf);
	}
		
	if ((gpflags & 16)){
		while (c=readByte()){
			//FILE COMMENT
		}
	}
	
	if ((gpflags & 2)){
		readByte();
		readByte();
	}
	
	DeflateLoop();
	
	crc = readByte();
	crc |= (readByte()<<8);
	crc |= (readByte()<<16);
	crc |= (readByte()<<24);
	
	size = readByte();
	size |= (readByte()<<8);
	size |= (readByte()<<16);
	size |= (readByte()<<24);
	
	nextFile();
	
};

};

/**
*  Base64 encoding / decoding
*  @see <a href="http://www.webtoolkit.info/">http://www.webtoolkit.info/</A>
*/
JXG.Util.Base64 = {

    // private property
    _keyStr : "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=",

    // public method for encoding
    encode : function (input) {
        var output = [],
            chr1, chr2, chr3, enc1, enc2, enc3, enc4,
            i = 0;

        input = JXG.Util.Base64._utf8_encode(input);

        while (i < input.length) {

            chr1 = input.charCodeAt(i++);
            chr2 = input.charCodeAt(i++);
            chr3 = input.charCodeAt(i++);

            enc1 = chr1 >> 2;
            enc2 = ((chr1 & 3) << 4) | (chr2 >> 4);
            enc3 = ((chr2 & 15) << 2) | (chr3 >> 6);
            enc4 = chr3 & 63;

            if (isNaN(chr2)) {
                enc3 = enc4 = 64;
            } else if (isNaN(chr3)) {
                enc4 = 64;
            }

            output.push([this._keyStr.charAt(enc1),
                         this._keyStr.charAt(enc2),
                         this._keyStr.charAt(enc3),
                         this._keyStr.charAt(enc4)].join(''));
        }

        return output.join('');
    },

    // public method for decoding
    decode : function (input, utf8) {
        var output = [],
            chr1, chr2, chr3,
            enc1, enc2, enc3, enc4,
            i = 0;

        input = input.replace(/[^A-Za-z0-9\+\/\=]/g, "");

        while (i < input.length) {

            enc1 = this._keyStr.indexOf(input.charAt(i++));
            enc2 = this._keyStr.indexOf(input.charAt(i++));
            enc3 = this._keyStr.indexOf(input.charAt(i++));
            enc4 = this._keyStr.indexOf(input.charAt(i++));

            chr1 = (enc1 << 2) | (enc2 >> 4);
            chr2 = ((enc2 & 15) << 4) | (enc3 >> 2);
            chr3 = ((enc3 & 3) << 6) | enc4;

            output.push(String.fromCharCode(chr1));

            if (enc3 != 64) {
                output.push(String.fromCharCode(chr2));
            }
            if (enc4 != 64) {
                output.push(String.fromCharCode(chr3));
            }
        }
        
        output = output.join(''); 
        
        if (utf8) {
            output = JXG.Util.Base64._utf8_decode(output);
        }
        return output;

    },

    // private method for UTF-8 encoding
    _utf8_encode : function (string) {
        string = string.replace(/\r\n/g,"\n");
        var utftext = "";

        for (var n = 0; n < string.length; n++) {

            var c = string.charCodeAt(n);

            if (c < 128) {
                utftext += String.fromCharCode(c);
            }
            else if((c > 127) && (c < 2048)) {
                utftext += String.fromCharCode((c >> 6) | 192);
                utftext += String.fromCharCode((c & 63) | 128);
            }
            else {
                utftext += String.fromCharCode((c >> 12) | 224);
                utftext += String.fromCharCode(((c >> 6) & 63) | 128);
                utftext += String.fromCharCode((c & 63) | 128);
            }

        }

        return utftext;
    },

    // private method for UTF-8 decoding
    _utf8_decode : function (utftext) {
        var string = [],
            i = 0,
            c = 0, c2 = 0, c3 = 0;

        while ( i < utftext.length ) {
            c = utftext.charCodeAt(i);
            if (c < 128) {
                string.push(String.fromCharCode(c));
                i++;
            }
            else if((c > 191) && (c < 224)) {
                c2 = utftext.charCodeAt(i+1);
                string.push(String.fromCharCode(((c & 31) << 6) | (c2 & 63)));
                i += 2;
            }
            else {
                c2 = utftext.charCodeAt(i+1);
                c3 = utftext.charCodeAt(i+2);
                string.push(String.fromCharCode(((c & 15) << 12) | ((c2 & 63) << 6) | (c3 & 63)));
                i += 3;
            }
        }
        return string.join('');
    },
    
    _destrip: function (stripped, wrap){
        var lines = [], lineno, i,
            destripped = [];
        
        if (wrap==null) 
            wrap = 76;
            
        stripped.replace(/ /g, "");
        lineno = stripped.length / wrap;
        for (i = 0; i < lineno; i++)
            lines[i]=stripped.substr(i * wrap, wrap);
        if (lineno != stripped.length / wrap)
            lines[lines.length]=stripped.substr(lineno * wrap, stripped.length-(lineno * wrap));
            
        for (i = 0; i < lines.length; i++)
            destripped.push(lines[i]);
        return destripped.join('\n');
    },
    
    decodeAsArray: function (input){
        var dec = this.decode(input),
            ar = [], i;
        for (i=0;i<dec.length;i++){
            ar[i]=dec.charCodeAt(i);
        }
        return ar;
    },
    
    decodeGEONExT : function (input) {
        return decodeAsArray(destrip(input),false);
    }
};

/**
 * @private
 */
JXG.Util.asciiCharCodeAt = function(str,i){
	var c = str.charCodeAt(i);
	if (c>255){
    	switch (c) {
			case 8364: c=128;
	    	break;
	    	case 8218: c=130;
	    	break;
	    	case 402: c=131;
	    	break;
	    	case 8222: c=132;
	    	break;
	    	case 8230: c=133;
	    	break;
	    	case 8224: c=134;
	    	break;
	    	case 8225: c=135;
	    	break;
	    	case 710: c=136;
	    	break;
	    	case 8240: c=137;
	    	break;
	    	case 352: c=138;
	    	break;
	    	case 8249: c=139;
	    	break;
	    	case 338: c=140;
	    	break;
	    	case 381: c=142;
	    	break;
	    	case 8216: c=145;
	    	break;
	    	case 8217: c=146;
	    	break;
	    	case 8220: c=147;
	    	break;
	    	case 8221: c=148;
	    	break;
	    	case 8226: c=149;
	    	break;
	    	case 8211: c=150;
	    	break;
	    	case 8212: c=151;
	    	break;
	    	case 732: c=152;
	    	break;
	    	case 8482: c=153;
	    	break;
	    	case 353: c=154;
	    	break;
	    	case 8250: c=155;
	    	break;
	    	case 339: c=156;
	    	break;
	    	case 382: c=158;
	    	break;
	    	case 376: c=159;
	    	break;
	    	default:
	    	break;
	    }
	}
	return c;
};

/**
 * Decoding string into utf-8
 * @param {String} string to decode
 * @return {String} utf8 decoded string
 */
JXG.Util.utf8Decode = function(utftext) {
  var string = [];
  var i = 0;
  var c = 0, c1 = 0, c2 = 0;

  while ( i < utftext.length ) {
    c = utftext.charCodeAt(i);

    if (c < 128) {
      string.push(String.fromCharCode(c));
      i++;
    } else if((c > 191) && (c < 224)) {
      c2 = utftext.charCodeAt(i+1);
      string.push(String.fromCharCode(((c & 31) << 6) | (c2 & 63)));
      i += 2;
    } else {
      c2 = utftext.charCodeAt(i+1);
      c3 = utftext.charCodeAt(i+2);
      string.push(String.fromCharCode(((c & 15) << 12) | ((c2 & 63) << 6) | (c3 & 63)));
      i += 3;
    }
  };
  return string.join('');
};
JXG.PsTricks = new function() {
    this.psTricksString = "";
};

JXG.PsTricks.convertBoardToPsTricks = function(board) {
    var p = new JXG.Coords(JXG.COORDS_BY_SCREEN, [board.canvasWidth, board.canvasHeight], board);
    var q = new JXG.Coords(JXG.COORDS_BY_SCREEN, [0, 0], board);
    this.psTricksString = '\\begin{pspicture*}('+q.usrCoords[1]+','+p.usrCoords[2]+')('+p.usrCoords[1]+','+q.usrCoords[2]+')\n';

    // Arcs (hier nur Sektoren)
    for(var el in board.objects) {
        var pEl = board.objects[el];
        if(pEl.type == JXG.OBJECT_TYPE_ARC) {
            if(pEl.visProp['visible']) {
                this.addSector(pEl);
            }
        }
    }    
    // Polygone
    for(var el in board.objects) {
        var pEl = board.objects[el];
        if(pEl.type == JXG.OBJECT_TYPE_POLYGON) {
            if(pEl.visProp['visible']) {
                this.addPolygon(pEl);
            }
        }
    }
    // Winkel
    for(var el in board.objects) {
        var pEl = board.objects[el];
        if(pEl.type == JXG.OBJECT_TYPE_ANGLE) {
            if(pEl.visProp['visible']) {
                this.addAngle(pEl);
            }
        }
    }
    // Kreise
    for(var el in board.objects) {
        var pEl = board.objects[el];
        if(pEl.type == JXG.OBJECT_TYPE_CIRCLE) {
            if(pEl.visProp['visible']) {
                this.addCircle(pEl);
            }
        }
    }
    // Arcs
    for(var el in board.objects) {
        var pEl = board.objects[el];
        if(pEl.type == JXG.OBJECT_TYPE_ARC) {
            if(pEl.visProp['visible']) {
                this.addArc(pEl);
            }
        }
    }
    // Linien
    for(var el in board.objects) {
        var pEl = board.objects[el];
        if(pEl.type == JXG.OBJECT_TYPE_LINE) {
            if(pEl.visProp['visible']) {
                this.addLine(pEl);
            }
        }
    }
    // Punkte
    for(var el in board.objects) {
        var pEl = board.objects[el];
        if(pEl.type == JXG.OBJECT_TYPE_POINT) {
            if(pEl.visProp['visible']) {
                this.addPoint(pEl);
            }
        }
    }    
    this.psTricksString += '\\end{pspicture*}';
};

JXG.PsTricks.givePsTricksToDiv = function(divId, board) {
    this.convertBoardToPsTricks(board);
    document.getElementById(divId).innerHTML = this.psTricksString;
};

JXG.PsTricks.addPoint = function(el) {
    this.psTricksString += "\\psdot";
    this.psTricksString += "[linecolor=" + this.parseColor(el.visProp['strokeColor']) + ",";
    this.psTricksString += "dotstyle=";
    if(el.visProp['face'] == 'cross') { // x
        this.psTricksString += "x, dotsize=";
        if(el.visProp['size'] == 2) {
            this.psTricksString += "2pt 2";
        }
        else if(el.visProp['size'] == 3) {
            this.psTricksString += "5pt 2";
        }
        else if(el.visProp['size'] >= 4) {
            this.psTricksString += "5pt 3";
        }        
    }
    else if(el.visProp['face'] == 'circle') { // circle
        this.psTricksString += "*, dotsize=";
        if(el.visProp['size'] == 1) {
            this.psTricksString += "2pt 2";
        }
        else if(el.visProp['size'] == 2) {
            this.psTricksString += "4pt 2";
        }
        else if(el.visProp['size'] == 3) {
            this.psTricksString += "6pt 2";
        }  
        else if(el.visProp['size'] >= 4) { // TODO
            this.psTricksString += "6pt 3";
        }          
    }
    else if(el.visProp['face'] == 'square') { // rectangle
        this.psTricksString += "square*, dotsize=";
        if(el.visProp['size'] == 2) {
            this.psTricksString += "2pt 2";
        }
        else if(el.visProp['size'] == 3) {
            this.psTricksString += "5pt 2";
        }
        else if(el.visProp['size'] >= 4) { // TODO
            this.psTricksString += "5pt 3";
        }           
    }
    else if(el.visProp['face'] == 'plus') { // +
        this.psTricksString += "+, dotsize=";
        if(el.visProp['size'] == 2) {
            this.psTricksString += "2pt 2";
        }
        else if(el.visProp['size'] == 3) {
            this.psTricksString += "5pt 2";
        }
        else if(el.visProp['size'] >= 4) { // TODO
            this.psTricksString += "5pt 3";
        }            
    }
    this.psTricksString += "]";
    this.psTricksString += "("+el.coords.usrCoords[1]+","+el.coords.usrCoords[2]+")\n";
    
    // Label
    this.psTricksString += "\\rput("+(el.coords.usrCoords[1]+15/ el.board.stretchY)+","+(el.coords.usrCoords[2]+15/ el.board.stretchY)+"){\\small $"+el.name+"$}\n";
};

JXG.PsTricks.addLine = function(el) {
    var screenCoords1 = new JXG.Coords(JXG.COORDS_BY_USER, el.point1.coords.usrCoords, el.board);
    var screenCoords2 = new JXG.Coords(JXG.COORDS_BY_USER, el.point2.coords.usrCoords, el.board);
    if(el.visProp['straightFirst'] || el.visProp['straightLast']) {
       el.board.renderer.calcStraight(el,screenCoords1,screenCoords2); 
    } 
    this.psTricksString += "\\psline";
    this.psTricksString += "[linecolor=" + this.parseColor(el.visProp['strokeColor']) + ", linewidth=" +el.visProp['strokeWidth']+"px";
    this.psTricksString += "]";
    if(el.visProp['firstArrow']) {
        if(el.visProp['lastArrow']) {
            this.psTricksString += "{<->}";
        }
        else {
            this.psTricksString += "{<-}";
        }
    }
    else {
        if(el.visProp['lastArrow']) {
            this.psTricksString += "{->}";
        }
    }
    this.psTricksString += "("+screenCoords1.usrCoords[1]+","+screenCoords1.usrCoords[2]+")("+screenCoords2.usrCoords[1]+","+screenCoords2.usrCoords[2]+")\n";
};

JXG.PsTricks.addCircle = function(el) {
    var radius = el.Radius();
    this.psTricksString += "\\pscircle";
    this.psTricksString += "[linecolor=" + this.parseColor(el.visProp['strokeColor']) +", linewidth=" +el.visProp['strokeWidth']+"px";
    if(el.visProp['fillColor'] != 'none' && el.visProp['fillOpacity'] != 0) {
        this.psTricksString += ", fillstyle=solid, fillcolor="+this.parseColor(el.visProp['fillColor'])+", opacity="+JXG.Math.round(el.visProp['fillOpacity'],5);
    }
    this.psTricksString += "]";
    this.psTricksString += "("+el.midpoint.coords.usrCoords[1]+","+el.midpoint.coords.usrCoords[2]+"){"+radius+"}\n";
};

JXG.PsTricks.addPolygon = function(el) {
    this.psTricksString += "\\pspolygon";
    this.psTricksString += "[linestyle=none, fillstyle=solid, fillcolor="+this.parseColor(el.visProp['fillColor'])+", opacity="+JXG.Math.round(el.visProp['fillOpacity'],5)+"]";
    for(var i=0; i < el.vertices.length; i++) {
        this.psTricksString += "("+el.vertices[i].coords.usrCoords[1]+","+el.vertices[i].coords.usrCoords[2]+")";
    }
    this.psTricksString += "\n";
};

JXG.PsTricks.addArc = function(el) {
    var radius = el.Radius();  
    var p = {};
    p.coords = new JXG.Coords(JXG.COORDS_BY_USER, 
                          [el.board.canvasWidth/(el.board.stretchY), el.midpoint.coords.usrCoords[2]],
                          el.board);
    var angle2 = JXG.Math.round(el.board.algebra.trueAngle(p, el.midpoint, el.point2),4);
    var angle1 = JXG.Math.round(el.board.algebra.trueAngle(p, el.midpoint, el.point3),4);
    
    this.psTricksString += "\\psarc";
    this.psTricksString += "[linecolor=" + this.parseColor(el.visProp['strokeColor']) + ", linewidth=" +el.visProp['strokeWidth']+"px";
    this.psTricksString += "]";
    if(el.visProp['lastArrow']) {
        if(el.visProp['firstArrow']) {
            this.psTricksString += "{<->}";
        }
        else {
            this.psTricksString += "{<-}";
        }
    }
    else {
        if(el.visProp['firstArrow']) {
            this.psTricksString += "{->}";
        }
    }    
    this.psTricksString += "("+el.midpoint.coords.usrCoords[1]+","+el.midpoint.coords.usrCoords[2]+"){"+radius+"}{"+angle2+"}{"+angle1+"}\n";
};

JXG.PsTricks.addSector = function(el) {
    var radius = el.Radius();  
    var p = {};
    p.coords = new JXG.Coords(JXG.COORDS_BY_USER, 
                          [el.board.canvasWidth/(el.board.stretchY), el.midpoint.coords.usrCoords[2]],
                          el.board);
    var angle2 = JXG.Math.round(el.board.algebra.trueAngle(p, el.midpoint, el.point2),4);
    var angle1 = JXG.Math.round(el.board.algebra.trueAngle(p, el.midpoint, el.point3),4);

    if(el.visProp['fillColor'] != 'none' && el.visProp['fillOpacity'] != 0) {
        this.psTricksString += "\\pswedge";
        this.psTricksString += "[linestyle=none, fillstyle=solid, fillcolor="+this.parseColor(el.visProp['fillColor'])+", opacity="+JXG.Math.round(el.visProp['fillOpacity'],5)+"]";
        this.psTricksString += "("+el.midpoint.coords.usrCoords[1]+","+el.midpoint.coords.usrCoords[2]+"){"+radius+"}{"+angle2+"}{"+angle1+"}\n";    
    }
};

JXG.PsTricks.addAngle = function(el) {
    var radius = el.radius;
    var p = {};
    p.coords = new JXG.Coords(JXG.COORDS_BY_USER, 
                          [el.board.canvasWidth/(el.board.stretchY), el.point2.coords.usrCoords[2]],
                          el.board);
    var angle2 = JXG.Math.round(el.board.algebra.trueAngle(p, el.point2, el.point1),4);
    var angle1 = JXG.Math.round(el.board.algebra.trueAngle(p, el.point2, el.point3),4);

    if(el.visProp['fillColor'] != 'none' && el.visProp['fillOpacity'] != 0) {
        this.psTricksString += "\\pswedge";
        this.psTricksString += "[linestyle=none, fillstyle=solid, fillcolor="+this.parseColor(el.visProp['fillColor'])+", opacity="+JXG.Math.round(el.visProp['fillOpacity'],5)+"]";
        this.psTricksString += "("+el.point2.coords.usrCoords[1]+","+el.point2.coords.usrCoords[2]+"){"+radius+"}{"+angle2+"}{"+angle1+"}\n";    
    }
    this.psTricksString += "\\psarc";
    this.psTricksString += "[linecolor=" + this.parseColor(el.visProp['strokeColor']) + ", linewidth=" +el.visProp['strokeWidth']+"px";
    this.psTricksString += "]"; 
    this.psTricksString += "("+el.point2.coords.usrCoords[1]+","+el.point2.coords.usrCoords[2]+"){"+radius+"}{"+angle2+"}{"+angle1+"}\n";
};

JXG.PsTricks.parseColor = function(color) {
    var arr = JXG.rgbParser(color);
    return "{[rgb]{"+arr[0]/255+","+arr[1]/255+","+arr[2]/255+"}}";
};

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph. If not, see <http://www.gnu.org/licenses/>.
*/

/** 
 * @fileoverview The JXG.Server is a wrapper for a smoother integration of server side calculations. on the
 * server side a python plugin system is used.
 */

/** TODO: Documentation */

/**
 * @namespace
 * JXG.Server namespace holding functions to load JXG server modules.
 */
JXG.Server = function(){};

/**
 * This is where all of a module's handlers are accessed from. If you're loading a module named JXGModule which
 * provides a handler called ImaHandler, then this handler can be called by invoking JXG.Server.modules.JXGModule.ImaHandler().
 * @namespace
 */
JXG.Server.modules = function(){};

/**
 * Stores all asynchronous calls to server which aren't finished yet.
 * @private
 */
JXG.Server.runningCalls = {};

/**
 * Handles errors, just a default implementation, can be overwritten by you, if you want to handle errors by yourself.
 * @param {object} data An object holding a field of type string named message handling the error described in the message string.
 */
JXG.Server.handleError = function(data) {
	alert('error occured, server says: ' + data.message);
};

/**
 * The main method of JXG.Server. Actually makes the calls to the server and parses the feedback.
 * @param {string} action Can be 'load' or 'exec'.
 * @param {function} callback Function pointer or anonymous function which takes as it's only argument an
 * object containing the data from the server. The fields of this object depend on the reply of the server
 * module. See the correspondings server module readme.
 * @param {object} data What is to be sent to the server.
 * @param {boolean} sync If the call should be synchronous or not.
 */
JXG.Server.callServer = function(action, callback, data, sync) {
	var fileurl, passdata, AJAX,
	params, id, dataJSONStr,
	k;

    if(typeof sync == 'undefined' || sync == null)
        sync = false;

	params = '';
	for(k in data) {
		params += '&' + escape(k) + '=' + escape(data[k]);
	}

	dataJSONStr = JXG.toJSON(data);

	// generate id
	do {
		id = action + Math.floor(Math.random()*4096);
	} while(typeof this.runningCalls[id] != 'undefined');

	// store information about the calls
	this.runningCalls[id] = { 'action': action };
	if(typeof data.module != 'undefined')
		this.runningCalls[id].module = data.module;

	fileurl = JXG.serverBase + 'JXGServer.py';
    passdata = 'action=' + escape(action) + '&id=' + id + '&dataJSON=' + escape(JXG.Util.Base64.encode(dataJSONStr));

	this.cbp = function(d) {
		var str, data,
		tmp, inject, paramlist, id,
		i, j;

		str = (new JXG.Util.Unzip(JXG.Util.Base64.decodeAsArray(d))).unzip();
		if(JXG.isArray(str) && str.length > 0)
			str = str[0][0];

        if(typeof str != 'string')
            return;

		data =  eval("(" + str + ")");

		if(data.type == 'error') {
			this.handleError(data);
		} else if (data.type == 'response') {
			id = data.id;

			// inject fields
			for(i=0; i<data.fields.length; i++) {
				tmp = data.fields[i];
				inject = tmp.namespace + ( typeof eval(tmp.namespace) == 'object' ? '.' : '.prototype.') + tmp.name + ' = ' + tmp.value;
				eval(inject);
			}

			// inject handlers
			for(i=0; i<data.handler.length; i++) {
				tmp = data.handler[i];
				paramlist = [];

				for(j=0; j<tmp.parameters.length; j++) {
					paramlist[j] = '"' + tmp.parameters[j] + '": ' + tmp.parameters[j];
				}
				// insert subnamespace named after module.
				inject = 'if(typeof JXG.Server.modules.' + this.runningCalls[id].module + ' == "undefined")' +
				'JXG.Server.modules.' + this.runningCalls[id].module + ' = {};';

				// insert callback method which fetches and uses the server's data for calculation in JavaScript
				inject += 'JXG.Server.modules.' + this.runningCalls[id].module + '.' + tmp.name + '_cb = ' + tmp.callback + ';';

				// insert handler as JXG.Server.modules.<module name>.<handler name>
				inject += 'JXG.Server.modules.' + this.runningCalls[id].module + '.' + tmp.name + ' = function (' + tmp.parameters.join(',') + ', __JXGSERVER_CB__) {' +
				'if(typeof __JXGSERVER_CB__ == "undefined") __JXGSERVER_CB__ = JXG.Server.modules.' + this.runningCalls[id].module + '.' + tmp.name + '_cb;' +
				'var __JXGSERVER_PAR__ = {' + paramlist.join(',') + ', "module": "' + this.runningCalls[id].module + '", "handler": "' + tmp.name + '" };' +
				'JXG.Server.callServer("exec", __JXGSERVER_CB__, __JXGSERVER_PAR__);' +
				'};';
				eval(inject);
			}

			delete this.runningCalls[id];

			// handle data
			callback(data.data);
		}
	};

	// bind cbp callback method to JXG.Server to get access to JXG.Server fields from within cpb
	this.cb = JXG.bind(this.cbp, this);

    // we're using our own XMLHttpRequest object in here because of a/sync and POST
    if (window.XMLHttpRequest) {
        AJAX = new XMLHttpRequest();
        AJAX.overrideMimeType('text/plain; charset=iso-8859-1');
    } else {                                  
        AJAX = new ActiveXObject("Microsoft.XMLHTTP");
    }
    if (AJAX) {
        // POST is required if data sent to server is too long for a url.
        // some browsers/http servers don't accept long urls.
        AJAX.open("POST", fileurl, !sync);
        AJAX.setRequestHeader("Content-type", "application/x-www-form-urlencoded");

        if(!sync) {
            // Define function to fetch data received from server
            // that function returning a function is required to make this.cb known to the function.
            AJAX.onreadystatechange = (function(cb){ return function () {
                switch(AJAX.readyState) {
                    // server is ready for take-off
                    case 4:
                        if(AJAX.status != 200)
                            alert("Fehler:" + AJAX.status);
                        else  // grab it and call the server callback to debase64, unzip, and parse the data
                            cb(AJAX.responseText);
                    break;
                    default:
                        return false;
                    break;
                }
            }})(this.cb);
        }

        // send the data
        AJAX.send(passdata);
        if(sync)
            this.cb(AJAX.responseText);
    } else {
        return false;
    }

//	JXG.FileReader.parseFileContent(fileurl, this.cb, 'raw', !sync);
};

/**
 * Callback for the default action 'load'.
 */
JXG.Server.loadModule_cb = function(data) {
	var i;
	for(i=0; i<data.length; i++)
		alert(data[i].name + ': ' + data[i].value);
};

/**
 * Loads a module from the server.
 * @param {string} module A string containing the module. Has to match the filename of the Python module on the server exactly including
 * lower and upper case letters without the file ending .py.
 */
JXG.Server.loadModule = function(module) {
	return JXG.Server.callServer('load', JXG.Server.loadModule_cb, {'module': module}, true);
};


/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph. If not, see <http://www.gnu.org/licenses/>.
*/

/** 
 * @fileoverview The JXG.DataSource is a helper class for data organization. Currently supported data sources are
 * javascript arrays and HTML tables.
 */

/* NOT YET DOCUMENTED. TODO! */

JXG.DataSource = function() {

    this.data = [];
    this.columnHeaders = [];
    this.rowHeaders = [];

    return this;
};

JXG.DataSource.prototype.loadFromArray = function(table, columnHeader, rowHeader) {
    var i, j,cell;

    if(typeof columnHeader == 'undefined')
        columnHeader = false;
    if(typeof rowHeader == 'undefined')
        rowHeader = false;

    if(JXG.isArray(columnHeader)) {
        this.columnHeader = columnHeader;
        columnHeader = false;
    }

    if(JXG.isArray(rowHeader)) {
        this.rowHeader = rowHeader;
        rowHeader = false;
    }

    this.data = [];
    if(columnHeader)
        this.columnHeader = [];
    if(rowHeader)
        this.rowHeader = [];

    if(typeof table != 'undefined') {
        // extract the data
        this.data = new Array(table.length);

        for(i=0; i<table.length; i++) {
            this.data[i] = new Array(table[i].length);
            for(j=0; j<table[i].length; j++) {
                cell = table[i][j];
                if('' + parseFloat(cell) == cell)
                    this.data[i][j] = parseFloat(cell);
                else if (cell != '-')
                    this.data[i][j] = cell;
                else
                    this.data[i][j] = NaN;
            }
        }
            
        if(columnHeader) {
            this.columnHeader = this.data[0].slice(1);
            this.data = this.data.slice(1);
        }

        if(rowHeader) {
            this.rowHeader = new Array();
            for(i=0; i<this.data.length; i++) {
                this.rowHeader.push(this.data[i][0]);
                this.data[i] = this.data[i].slice(1);
            }
        }
    }

    return this;
};

JXG.DataSource.prototype.loadFromTable = function(table, columnHeader, rowHeader) {
    var row, i, j, col, cell, name;

    if(typeof columnHeader == 'undefined')
        columnHeader = false;
    if(typeof rowHeader == 'undefined')
        rowHeader = false;

    if(JXG.isArray(columnHeader)) {
        this.columnHeader = columnHeader;
        columnHeader = false;
    }

    if(JXG.isArray(rowHeader)) {
        this.rowHeader = rowHeader;
        rowHeader = false;
    }

    this.data = [];
    if(columnHeader)
        this.columnHeader = [];
    if(rowHeader)
        this.rowHeader = [];

    table = document.getElementById(table);
    if(typeof table != 'undefined') {
        // extract the data
        row = table.getElementsByTagName('tr');
        this.data = new Array(row.length);

        for(i=0; i<row.length; i++) {
            col = row[i].getElementsByTagName('td');
            this.data[i] = new Array(col.length);
            for(j=0; j<col.length; j++) {
                cell = col[j].innerHTML;
                if('' + parseFloat(cell) == cell)
                    this.data[i][j] = parseFloat(cell);
                else if (cell != '-')
                    this.data[i][j] = cell;
                else
                    this.data[i][j] = NaN;
            }
        }
            
        if(columnHeader) {
            this.columnHeader = this.data[0].slice(1);
            this.data = this.data.slice(1);
        }

        if(rowHeader) {
            this.rowHeader = new Array();
            for(i=0; i<this.data.length; i++) {
                this.rowHeader.push(this.data[i][0]);
                this.data[i] = this.data[i].slice(1);
            }
        }
    }

    return this;
};

JXG.DataSource.prototype.addColumn = function(name, pos, data) {
    // todo
};

JXG.DataSource.prototype.addRow = function(name, pos, data) {
    // todo
};

JXG.DataSource.prototype.getColumn = function(col) {
    var result = new Array(this.data.length), i;

    // get column index if column is given as column header title
    if(typeof col == 'string') {
        for(i=0; i<this.columnHeader.length; i++) {
            if(col == this.columnHeader[i]) {
                col = i;
                break;
            }
        }
    }

    // build column array
    for(i=0; i<this.data.length; i++) {
        if(this.data[i].length > col)
            result[i] = this.data[i][col];
    }

    return result;
};

JXG.DataSource.prototype.getRow = function(row) {
    var result, i;

    // get column index if column is given as column header title
    if(typeof row == 'string') {
        for(i=0; i<this.rowHeader.length; i++) {
            if(row == this.rowHeader[i]) {
                row = i;
                break;
            }
        }
    }

    // allocate memory for result array
    result = new Array(this.data[row].length);

    // build column array. result = this.data[row] is a flat copy and will
    // destroy our local data copy, that's why we're copying it element wise.
    for(i=0; i<this.data[row].length; i++) {
        result[i] = this.data[row][i];
    }

    return result;
};

/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

JXG.SVGRenderer = function(container) {
    var i;
    this.constructor();

    this.svgRoot = null;
    this.suspendHandle = null;
    
    this.svgNamespace = 'http://www.w3.org/2000/svg';
    this.xlinkNamespace ='http://www.w3.org/1999/xlink';

    this.container = container;
    this.container.style.MozUserSelect = 'none';

    this.container.style.overflow = 'hidden';
    if (this.container.style.position=='') {
        this.container.style.position = 'relative';
    }
    
    this.svgRoot = this.container.ownerDocument.createElementNS(this.svgNamespace, "svg");
    this.container.appendChild(this.svgRoot);

    this.defs = this.container.ownerDocument.createElementNS(this.svgNamespace,'defs');
    this.svgRoot.appendChild(this.defs);
    this.filter = this.container.ownerDocument.createElementNS(this.svgNamespace,'filter');
    this.filter.setAttributeNS(null, 'id', 'f1');
    this.filter.setAttributeNS(null, 'width', '300%');
    this.filter.setAttributeNS(null, 'height', '300%');
    this.feOffset = this.container.ownerDocument.createElementNS(this.svgNamespace,'feOffset');
    this.feOffset.setAttributeNS(null, 'result', 'offOut');
    this.feOffset.setAttributeNS(null, 'in', 'SourceAlpha');
    this.feOffset.setAttributeNS(null, 'dx', '5');
    this.feOffset.setAttributeNS(null, 'dy', '5');
    this.filter.appendChild(this.feOffset);
    this.feGaussianBlur = this.container.ownerDocument.createElementNS(this.svgNamespace,'feGaussianBlur');
    this.feGaussianBlur.setAttributeNS(null, 'result', 'blurOut');
    this.feGaussianBlur.setAttributeNS(null, 'in', 'offOut');
    this.feGaussianBlur.setAttributeNS(null, 'stdDeviation', '3');
    this.filter.appendChild(this.feGaussianBlur);
    this.feBlend = this.container.ownerDocument.createElementNS(this.svgNamespace,'feBlend');
    this.feBlend.setAttributeNS(null, 'in', 'SourceGraphic');
    this.feBlend.setAttributeNS(null, 'in2', 'blurOut');
    this.feBlend.setAttributeNS(null, 'mode', 'normal');
    this.filter.appendChild(this.feBlend);
    this.defs.appendChild(this.filter);    
    
    // um eine passende Reihenfolge herzustellen
    /*
    this.images = this.container.ownerDocument.createElementNS(this.svgNamespace,'g');
    this.svgRoot.appendChild(this.images);
    this.grid = this.container.ownerDocument.createElementNS(this.svgNamespace,'g');
    this.svgRoot.appendChild(this.grid);
    this.angles = this.container.ownerDocument.createElementNS(this.svgNamespace,'g');
    this.svgRoot.appendChild(this.angles);    
    this.sectors = this.container.ownerDocument.createElementNS(this.svgNamespace,'g');
    this.svgRoot.appendChild(this.sectors);
    this.polygone = this.container.ownerDocument.createElementNS(this.svgNamespace,'g');
    this.svgRoot.appendChild(this.polygone);
    this.curves = this.container.ownerDocument.createElementNS(this.svgNamespace,'g');
    this.svgRoot.appendChild(this.curves);
    this.circles = this.container.ownerDocument.createElementNS(this.svgNamespace,'g');
    this.svgRoot.appendChild(this.circles);
    this.lines = this.container.ownerDocument.createElementNS(this.svgNamespace,'g');
    this.svgRoot.appendChild(this.lines);
    this.arcs = this.container.ownerDocument.createElementNS(this.svgNamespace,'g');
    this.svgRoot.appendChild(this.arcs);
    this.points = this.container.ownerDocument.createElementNS(this.svgNamespace,'g');
    this.svgRoot.appendChild(this.points);
    */
    /* 
    * 10 Layers. highest number = highest visibility
    */
    this.layer = [];
    for (i=0;i<JXG.Options.layer.numlayers;i++) {
        this.layer[i] = this.container.ownerDocument.createElementNS(this.svgNamespace,'g');
        this.svgRoot.appendChild(this.layer[i]);
    }
    
    // um Dashes zu realisieren
    this.dashArray = ['2, 2', '5, 5', '10, 10', '20, 20', '20, 10, 10, 10', '20, 5, 10, 5'];
};

JXG.SVGRenderer.prototype = new JXG.AbstractRenderer;

JXG.SVGRenderer.prototype.setShadow = function(el) {
    if (el.visPropOld['shadow']==el.visProp['shadow']) {
        return;
    }
    if(el.rendNode != null) {
        if(el.visProp['shadow']) {
            el.rendNode.setAttributeNS(null,'filter','url(#f1)');
        }
        else {
            el.rendNode.removeAttributeNS(null,'filter');
        }    
    }
    el.visPropOld['shadow']=el.visProp['shadow'];
}

JXG.SVGRenderer.prototype.setGradient = function(el) {
    var fillNode = el.rendNode, col, op;
    
    if(el.type == JXG.OBJECT_TYPE_ARC || el.type == JXG.OBJECT_TYPE_ANGLE) {
        fillNode = el.rendNode2;
    } 
    if (typeof el.visProp['fillOpacity']=='function') {
        op = el.visProp['fillOpacity']();
    } else {
        op = el.visProp['fillOpacity'];
    }
    op = (op>0)?op:0;
    if (typeof el.visProp['fillColor']=='function') {
        col = el.visProp['fillColor']();
    } else {
        col = el.visProp['fillColor'];
    }

    if(el.visProp['gradient'] == 'linear') {
        var node = this.createPrimitive('linearGradient',el.id+'_gradient');
        var x1 = '0%'; // TODO: get x1,x2,y1,y2 from el.visProp['angle']
        var x2 = '100%';
        var y1 = '0%';
        var y2 = '0%'; //means 270 degrees

        node.setAttributeNS(null,'x1',x1);
        node.setAttributeNS(null,'x2',x2);
        node.setAttributeNS(null,'y1',y1);
        node.setAttributeNS(null,'y2',y2);
        var node2 = this.createPrimitive('stop',el.id+'_gradient1');
        node2.setAttributeNS(null,'offset','0%');
        node2.setAttributeNS(null,'style','stop-color:'+col+';stop-opacity:'+op);     
        var node3 = this.createPrimitive('stop',el.id+'_gradient2');
        node3.setAttributeNS(null,'offset','100%');
        node3.setAttributeNS(null,'style','stop-color:'+el.visProp['gradientSecondColor']+';stop-opacity:'+el.visProp['gradientSecondOpacity']);
        node.appendChild(node2);
        node.appendChild(node3);     
        this.defs.appendChild(node);
        fillNode.setAttributeNS(null, 'style', 'fill:url(#'+el.id+'_gradient)');      
        el.gradNode1 = node2;
        el.gradNode2 = node3;
    }
    else if (el.visProp['gradient'] == 'radial') {
        var node = this.createPrimitive('radialGradient',el.id+'_gradient');

        node.setAttributeNS(null, 'cx', '50%')
        node.setAttributeNS(null, 'cy', '50%')
        node.setAttributeNS(null, 'r', '50%')
        node.setAttributeNS(null, 'fx', el.visProp['gradientPositionX']*100+'%')
        node.setAttributeNS(null, 'fy', el.visProp['gradientPositionY']*100+'%')

        var node2 = this.createPrimitive('stop',el.id+'_gradient1');
        node2.setAttributeNS(null,'offset','0%');
        node2.setAttributeNS(null,'style','stop-color:'+el.visProp['gradientSecondColor']+';stop-opacity:'+el.visProp['gradientSecondOpacity']);
        var node3 = this.createPrimitive('stop',el.id+'_gradient2');
        node3.setAttributeNS(null,'offset','100%');
        node3.setAttributeNS(null,'style','stop-color:'+col+';stop-opacity:'+op);         

        node.appendChild(node2);
        node.appendChild(node3);     
        this.defs.appendChild(node);
        fillNode.setAttributeNS(null, 'style', 'fill:url(#'+el.id+'_gradient)'); 
        el.gradNode1 = node2;
        el.gradNode2 = node3;
    }
    else {
        fillNode.removeAttributeNS(null,'style');
    }
};

JXG.SVGRenderer.prototype.updateGradient = function(el) {
    var node2 = el.gradNode1, 
        node3 = el.gradNode2, 
        col, op;

    if (node2==null || node3==0) {
        return;
    }
    if (typeof el.visProp['fillOpacity']=='function') {
        op = el.visProp['fillOpacity']();
    } else {
        op = el.visProp['fillOpacity'];
    }
    op = (op>0)?op:0;
    if (typeof el.visProp['fillColor']=='function') {
        col = el.visProp['fillColor']();
    } else {
        col = el.visProp['fillColor'];
    }
    
    if(el.visProp['gradient'] == 'linear') {
        node2.setAttributeNS(null,'style','stop-color:'+col+';stop-opacity:'+op);     
        node3.setAttributeNS(null,'style','stop-color:'+el.visProp['gradientSecondColor']+';stop-opacity:'+el.visProp['gradientSecondOpacity']);
    } else if (el.visProp['gradient'] == 'radial') {
        node2.setAttributeNS(null,'style','stop-color:'+el.visProp['gradientSecondColor']+';stop-opacity:'+el.visProp['gradientSecondOpacity']);
        node3.setAttributeNS(null,'style','stop-color:'+col+';stop-opacity:'+op);         
    }
}; 

JXG.SVGRenderer.prototype.displayCopyright = function(str,fontsize) {
    var node = this.createPrimitive('text','licenseText'),
        t;
    node.setAttributeNS(null,'x','20');
    node.setAttributeNS(null,'y',2+fontsize);
    node.setAttributeNS(null, "style", "font-family:Arial,Helvetica,sans-serif; font-size:"+fontsize+"px; fill:#356AA0;  opacity:0.3;");
    t = document.createTextNode(str);
    node.appendChild(t);
    this.appendChildPrimitive(node,0);
};

JXG.SVGRenderer.prototype.drawInternalText = function(el) {
    var node = this.createPrimitive('text',el.id);
    node.setAttributeNS(null, "class", "JXGtext");
    el.rendNodeText = document.createTextNode('');
    node.appendChild(el.rendNodeText);
    this.appendChildPrimitive(node,9);
    return node;
};

JXG.SVGRenderer.prototype.updateInternalText = function(/** JXG.Text */ el) { 
    el.rendNode.setAttributeNS(null, 'x', (el.coords.scrCoords[1])+'px'); 
    el.rendNode.setAttributeNS(null, 'y', (el.coords.scrCoords[2] - this.vOffsetText)+'px'); 
    el.updateText();
    if (el.htmlStr!= el.plaintextStr) {
        el.rendNodeText.data = el.plaintextStr;
        el.htmlStr = el.plaintextStr;
    }
};

JXG.SVGRenderer.prototype.drawTicks = function(axis) {
    var node = this.createPrimitive('path', axis.id);
    //node.setAttributeNS(null, 'shape-rendering', 'crispEdges');
    this.appendChildPrimitive(node,axis.layer);
    this.appendNodesToElement(axis,'path'); 
};

JXG.SVGRenderer.prototype.updateTicks = function(axis,dxMaj,dyMaj,dxMin,dyMin) {
    var tickStr = "",
        i, c, node, 
        len = axis.ticks.length;
        
    for (i=0; i<len; i++) {
        c = axis.ticks[i].scrCoords;
        if (axis.ticks[i].major) {
            if (axis.labels[i].visProp['visible']) this.drawText(axis.labels[i]);
            tickStr += "M " + (c[1]+dxMaj) + " " + (c[2]-dyMaj) + " L " + (c[1]-dxMaj) + " " + (c[2]+dyMaj) + " ";
        }
        else
            tickStr += "M " + (c[1]+dxMin) + " " + (c[2]-dyMin) + " L " + (c[1]-dxMin) + " " + (c[2]+dyMin) + " ";
    }
    
    node = document.getElementById(axis.id);
    if(node == null) {
        node = this.createPrimitive('path', axis.id);
        //node.setAttributeNS(null, 'shape-rendering', 'crispEdges');
        this.appendChildPrimitive(node,axis.layer);
        this.appendNodesToElement(axis,'path');
    }
    node.setAttributeNS(null, 'stroke', axis.visProp['strokeColor']);    
    node.setAttributeNS(null, 'stroke-opacity', axis.visProp['strokeOpacity']);
    node.setAttributeNS(null, 'stroke-width', axis.visProp['strokeWidth']);
    this.updatePathPrimitive(node, tickStr, axis.board);
};

JXG.SVGRenderer.prototype.drawArc = function(el) {  
    var node = this.createPrimitive('path',el.id),
        radius, angle, circle, point3,
        pathString, pathString2, node2, node4;
        
    el.rendNode = node;

    JXG.clearVisPropOld(el);
    
    radius = el.Radius();  
    angle = el.board.algebra.trueAngle(el.point2, el.midpoint, el.point3);
    circle = {}; // um projectToCircle benutzen zu koennen...
    circle.midpoint = el.midpoint;
    circle.Radius = function() {
        return radius;
    };
    //-------------------
    // deprecated
    circle.getRadius = function() {
        return radius;
    };
    //-------------------
    point3 = el.board.algebra.projectPointToCircle(el.point3,circle);

    pathString = 'M '+ el.point2.coords.scrCoords[1] +' '+ el.point2.coords.scrCoords[2] +' A '; // Startpunkt
    pathString += Math.round(radius * el.board.stretchX) + ' ' + Math.round(radius * el.board.stretchY) + ' 0 '; // Radien
    // largeArc
    if(angle >= 180) {
        pathString += '1 ';
    }
    else {
        pathString += '0 ';
    }
    // sweepFlag
    pathString += '0 ';
    pathString += point3.scrCoords[1] + ' ' + point3.scrCoords[2]; // Endpunkt
    
    this.updatePathPrimitive(node,pathString,el.board);
    if (el.visProp['strokeColor']!=null) {node.setAttributeNS(null, 'stroke', el.visProp['strokeColor']);}
    if (el.visProp['strokeOpacity']!=null) {node.setAttributeNS(null, 'stroke-opacity', el.visProp['strokeOpacity']);}
    if (el.visProp['strokeWidth']!=null) {node.setAttributeNS(null, 'stroke-width', el.visProp['strokeWidth']);}
    node.setAttributeNS(null, 'fill', 'none');
    this.setDashStyle(el,el.visProp);
    
    this.setShadow(el);

    if(el.visProp['firstArrow']) {
        node2 = this.createArrowHead(el,'Start');
        this.defs.appendChild(node2);
        el.rendNodeTriangleStart = node2;
        node.setAttributeNS(null, 'marker-end', 'url(#'+el.id+'TriangleStart)');
    }
    if(el.visProp['lastArrow']) {
        node2 = this.createArrowHead(el,'End');
        this.defs.appendChild(node2);
        el.rendNodeTriangleEnd = node2;
        node.setAttributeNS(null, 'marker-start', 'url(#'+el.id+'TriangleEnd)');
    }      
    
    // Fuellflaeche
    node4 = this.createPrimitive('path',el.id+'sector');
    el.rendNode2 = node4;
    
    pathString2 = 'M ' + el.midpoint.coords.scrCoords[1] + " " + el.midpoint.coords.scrCoords[2];
    pathString2 += ' L '+ el.point2.coords.scrCoords[1] +' '+ el.point2.coords.scrCoords[2] +' A '; // Startpunkt
    pathString2 += Math.round(radius * el.board.stretchX) + ' ' + Math.round(radius * el.board.stretchY) + ' 0 '; // Radien
    // largeArc
    if(angle >= 180) {
        pathString2 += '1 ';
    }
    else {
        pathString2 += '0 ';
    }
    // sweepFlag
    pathString2 += '0 ';
    pathString2 += point3.scrCoords[1] + ' ' + point3.scrCoords[2];
    pathString2 += ' L ' + el.midpoint.coords.scrCoords[1] + " " + el.midpoint.coords.scrCoords[2]    + ' z'; // Endpunkt
    
    this.updatePathPrimitive(node4,pathString2,el.board);
    if (el.visProp['fillColor']!=null) {node4.setAttributeNS(null, 'fill', el.visProp['fillColor']);}
    if (el.visProp['fillOpacity']!=null) {node4.setAttributeNS(null, 'fill-opacity', el.visProp['fillOpacity']);}     
    node4.setAttributeNS(null, 'stroke', 'none');
    this.setGradient(el);
    
    this.appendChildPrimitive(node,el.layer); 
    this.appendChildPrimitive(node4,2); // hard coded layer
    /*
    this.arcs.appendChild(node);
    this.sectors.appendChild(node4);
    */

    if (el.visProp['draft']) {
        this.setDraft(el);
    }
    if(!el.visProp['visible']) {
        el.hideElement();
    }
};

/**
 * Updates properties of an arc that already exists.
 * @param {JXG.Arc} arc Reference to an arc object, that has to be updated.
 * @see JXG.Arc
 * @see #drawArc
 */
JXG.SVGRenderer.prototype.updateArc = function(el) { 
    // brutaler Fix der update-Methode...
    var node;
   
    this.remove(el.rendNode);
    this.remove(el.rendNode2);     
    node = el.rendNodeTriangleStart;
    if (node != null) {
        this.remove(node);
    }
    node = el.rendNodeTriangleEnd;
    if (node != null) {
        this.remove(node);
    }    
    this.drawArc(el);
    return;
};

JXG.SVGRenderer.prototype.drawAngle = function(el) {
    var angle = el.board.algebra.trueAngle(el.point1, el.point2, el.point3),
        circle, projectedP1, projectedP3,
        node, node2, pathString;
    JXG.clearVisPropOld(el);

    circle = {};  // um projectToCircle benutzen zu koennen...
    circle.midpoint = el.point2;
    circle.Radius = function() {
        return el.radius;
    };
    //-------------------
    // deprecated
    circle.getRadius = function() {
        return el.radius;
    };
    //-------------------
    projectedP1 = el.board.algebra.projectPointToCircle(el.point1,circle);
    projectedP3 = el.board.algebra.projectPointToCircle(el.point3,circle);

    node = this.createPrimitive('path',el.id+'_1');
    pathString = 'M ' + el.point2.coords.scrCoords[1] + " " + el.point2.coords.scrCoords[2];
    pathString += ' L '+ projectedP1.scrCoords[1] +' '+ projectedP1.scrCoords[2] +' A '; // Startpunkt
    pathString += Math.round(el.radius * el.board.stretchX) + ' ' + Math.round(el.radius * el.board.stretchY) + ' 0 '; // Radien
    // largeArc
    if (angle >= 180) {
        pathString += '1 ';
    }
    else {
        pathString += '0 ';
    }
    // sweepFlag
    pathString += '0 ';
    pathString += projectedP3.scrCoords[1] + ' ' + projectedP3.scrCoords[2];
    pathString += ' L ' + el.point2.coords.scrCoords[1] + " " + el.point2.coords.scrCoords[2]    + ' z'; // Endpunkt

    node.setAttributeNS(null, 'd', pathString);    
    
    node.setAttributeNS(null, 'fill', el.visProp['fillColor']);
    node.setAttributeNS(null, 'fill-opacity', el.visProp['fillOpacity']);    
    node.setAttributeNS(null, 'stroke', 'none');    
   
    node2 = this.createPrimitive('path',el.id+'_2');
    pathString = 'M '+  projectedP1.scrCoords[1] +' '+  projectedP1.scrCoords[2] +' A '; // Startpunkt
    pathString += Math.round(el.radius * el.board.stretchX) + ' ' + Math.round(el.radius * el.board.stretchY) + ' 0 '; // Radien
    // largeArc
    if (angle >= 180) {
        pathString += '1 ';
    }
    else {
        pathString += '0 ';
    }
    // sweepFlag
    pathString += '0 ';
    pathString += projectedP3.scrCoords[1] + ' ' + projectedP3.scrCoords[2]; // Endpunkt    

    node2.setAttributeNS(null, 'd', pathString);
    node2.setAttributeNS(null, 'id', el.id+'_2');
    node2.setAttributeNS(null, 'fill', 'none');    
    node2.setAttributeNS(null, 'stroke', el.visProp['strokeColor']);    
    node2.setAttributeNS(null, 'stroke-opacity', el.visProp['strokeOpacity']);

    this.appendChildPrimitive(node,el.layer);
    el.rendNode = node;
    this.setShadow(el);    
    this.appendChildPrimitive(node2,2);  // hard coded layer
    el.rendNode2 = node2;

    this.setObjectStrokeWidth(el,el.visProp['strokeWidth']);
};

JXG.SVGRenderer.prototype.updateAngle = function(el) {
    /* erstmal nur der brutale Weg... */
    this.remove(el.rendNode);
    this.remove(el.rendNode2);    
    this.drawAngle(el);
    if (!el.visProp['visible']) {
        el.hideElement();
    }
    return;
};

JXG.SVGRenderer.prototype.drawImage = function(el) {
    var url = el.url, //'data:image/png;base64,' + el.imageBase64String,    
        node = this.createPrimitive('image',el.id);

    node.setAttributeNS(this.xlinkNamespace, 'xlink:href', url);
    node.setAttributeNS(null, 'preserveAspectRatio', 'none');
    this.appendChildPrimitive(node,el.layer);
    el.rendNode = node;
    this.updateImage(el);
};

JXG.SVGRenderer.prototype.transformImage = function(el,t) {
    var node = el.rendNode,
        str = node.getAttributeNS(null, 'transform');
        
    str += ' ' + this.joinTransforms(el,t);
    node.setAttributeNS(null, 'transform', str);
};

JXG.SVGRenderer.prototype.joinTransforms = function(el,t) {
    var str = '', i, s,
        len = t.length;
        
    for (i=0;i<len;i++) {
        s = t[i].matrix[1][1]+','+t[i].matrix[2][1]+','+t[i].matrix[1][2]+','+t[i].matrix[2][2]+','+t[i].matrix[1][0]+','+t[i].matrix[2][0];
        str += 'matrix('+s+') ';
    }
    return str;
};
  
JXG.SVGRenderer.prototype.transformImageParent = function(el,m) {
    var s, str;
    if (m!=null) {
        s = m[1][1]+','+m[2][1]+','+m[1][2]+','+m[2][2]+','+m[1][0]+','+m[2][0];
        str = 'matrix('+s+')';
    } else {
        str = '';
    }
    el.rendNode.setAttributeNS(null, 'transform', str);
};
  
/*
JXG.SVGRenderer.prototype.removeGrid = function(board) { 
    var c = this.layer[board.options.layer['grid']];
    board.hasGrid = false;
    while (c.childNodes.length>0) {
        c.removeChild(c.firstChild);
    }
};
*/
 
JXG.SVGRenderer.prototype.setObjectStrokeColor = function(el, color, opacity) {
    var c = this.eval(color), 
        o = this.eval(opacity), 
        node;

    o = (o>0)?o:0;

    if (el.visPropOld['strokeColor']==c && el.visPropOld['strokeOpacity']==o) {
        return;
    }
    node = el.rendNode;
    if(el.type == JXG.OBJECT_TYPE_TEXT) {
        node.style.color = c; // Schriftfarbe
    }
    else {
        node.setAttributeNS(null, 'stroke', c);
        node.setAttributeNS(null, 'stroke-opacity', o);          
    }
    if(el.type == JXG.OBJECT_TYPE_ARROW) {
         el.rendNodeTriangle.setAttributeNS(null, 'stroke', c);
         el.rendNodeTriangle.setAttributeNS(null, 'stroke-opacity', o);
         el.rendNodeTriangle.setAttributeNS(null, 'fill', c);
         el.rendNodeTriangle.setAttributeNS(null, 'fill-opacity', o);             
    }
    if(el.type == JXG.OBJECT_TYPE_ARC) {
        if(el.visProp['firstArrow']) {
            el.rendNodeTriangleStart.setAttributeNS(null, 'stroke', c);
            el.rendNodeTriangleStart.setAttributeNS(null, 'stroke-opacity', o);                
            el.rendNodeTriangleStart.setAttributeNS(null, 'fill', c);
            el.rendNodeTriangleStart.setAttributeNS(null, 'fill-opacity', o);                    
        }
        if(el.visProp['lastArrow']) {
            el.rendNodeTriangleEnd.setAttributeNS(null, 'stroke', c);
            el.rendNodeTriangleEnd.setAttributeNS(null, 'stroke-opacity', o);                
            el.rendNodeTriangleEnd.setAttributeNS(null, 'fill', c);
            el.rendNodeTriangleEnd.setAttributeNS(null, 'fill-opacity', o);    
        }                
    }     
    else if(el.type == JXG.OBJECT_TYPE_LINE) {
        if(el.visProp['firstArrow']) {
            el.rendNodeTriangleStart.setAttributeNS(null, 'stroke', c);
            el.rendNodeTriangleStart.setAttributeNS(null, 'stroke-opacity', o);                
            el.rendNodeTriangleStart.setAttributeNS(null, 'fill', c);
            el.rendNodeTriangleStart.setAttributeNS(null, 'fill-opacity', o);                    
        }
        if(el.visProp['lastArrow']) {
            el.rendNodeTriangleEnd.setAttributeNS(null, 'stroke', c);
            el.rendNodeTriangleEnd.setAttributeNS(null, 'stroke-opacity', o);                
            el.rendNodeTriangleEnd.setAttributeNS(null, 'fill', c);
            el.rendNodeTriangleEnd.setAttributeNS(null, 'fill-opacity', o);    
        }                
    }
    el.visPropOld['strokeColor'] = c;
    el.visPropOld['strokeOpacity'] = o;
};

JXG.SVGRenderer.prototype.setObjectFillColor = function(el, color, opacity) {
    var c = this.eval(color), 
        o = this.eval(opacity);

    o = (o>0)?o:0;

    if (el.visPropOld['fillColor']==c && el.visPropOld['fillOpacity']==o) {
        return;
    }
    if(el.type == JXG.OBJECT_TYPE_ARC || el.type == JXG.OBJECT_TYPE_ANGLE) {
        node = el.rendNode2;
        node.setAttributeNS(null, 'fill', c);
        node.setAttributeNS(null, 'fill-opacity', o);        
    }    
    else {
        node = el.rendNode;
        node.setAttributeNS(null, 'fill', c);           
        node.setAttributeNS(null, 'fill-opacity', o);                   
    }
    
    if (el.visProp['gradient']!=null) {
        this.updateGradient(el);
    }
    el.visPropOld['fillColor'] = c;
    el.visPropOld['fillOpacity'] = o;
} ;

/**
 * Sets an elements stroke width.
 * @param {Object} el Reference to the geometry element.
 * @param {int} width The new stroke width to be assigned to the element.
 */
JXG.SVGRenderer.prototype.setObjectStrokeWidth = function(el, width) {
    var w = this.eval(width), 
        node;
    //w = (w>0)?w:0;
    try {
        if (el.visPropOld['strokeWidth']==w) {
            return;
        }
    } catch (e){
        //alert(el.id);
    }
    
    if(el.elementClass != JXG.OBJECT_CLASS_POINT) {
        if(el.type == JXG.OBJECT_TYPE_ANGLE) {
            node = el.rendNode2;
        }
        else {
            node = el.rendNode;
        }
        this.setPropertyPrimitive(node,'stroked', 'true');
        if (w!=null) { 
            this.setPropertyPrimitive(node,'stroke-width',w);    
        }
    }
    else {
        node = el.rendNode;
        this.setPropertyPrimitive(node,'stroked', 'true');
        if (w!=null) { 
            this.setPropertyPrimitive(node,'stroke-width',w); 
        }
    }
    el.visPropOld['strokeWidth'] = w;
};

JXG.SVGRenderer.prototype.hide = function(el) {
    var node;
    if (el==null) return;
    if(el.type == JXG.OBJECT_TYPE_ARC) {
        node = el.rendNode;
        node.setAttributeNS(null, 'display', 'none');
        node.style.visibility = "hidden"; 
        node = el.rendNode2;
        node.setAttributeNS(null, 'display', 'none');
        node.style.visibility = "hidden";         
    }
    else if(el.type == JXG.OBJECT_TYPE_ANGLE) {
        node = el.rendNode;
        node.setAttributeNS(null, 'display', 'none');
        node.style.visibility = "hidden"; 
        node = el.rendNode2;
        node.setAttributeNS(null, 'display', 'none');
        node.style.visibility = "hidden";         
    }   
    else {
        node = el.rendNode;
        node.setAttributeNS(null, 'display', 'none');
        node.style.visibility = "hidden";     
    }
};

JXG.SVGRenderer.prototype.show = function(el) {
    var node;
    if(el.type == JXG.OBJECT_TYPE_ARC) {
        node = el.rendNode;
        node.setAttributeNS(null, 'display', 'inline');
        node.style.visibility = "inherit"; 
        node = el.rendNode2;
        node.setAttributeNS(null, 'display', 'inline');
        node.style.visibility = "inherit";     
    }
    else if(el.type == JXG.OBJECT_TYPE_ANGLE) {
        node = el.rendNode;
        node.setAttributeNS(null, 'display', 'inline');
        node.style.visibility = "inherit"; 
        node = el.rendNode2;
        node.setAttributeNS(null, 'display', 'inline');
        node.style.visibility = "inherit";         
    }    
    else {
        node = el.rendNode;
        node.setAttributeNS(null, 'display', 'inline');
        node.style.visibility = "inherit"; 
    }
};

JXG.SVGRenderer.prototype.remove = function(shape) {
    if(shape!=null && shape.parentNode != null)
        shape.parentNode.removeChild(shape);
};

JXG.SVGRenderer.prototype.suspendRedraw = function() {
    // It seems to be important for the Linux version of firefox
    if (true) { this.suspendHandle = this.svgRoot.suspendRedraw(10000); }
};

JXG.SVGRenderer.prototype.unsuspendRedraw = function() {
    if (true) { 
        this.svgRoot.unsuspendRedraw(this.suspendHandle);
        this.svgRoot.forceRedraw();
    }
};

JXG.SVGRenderer.prototype.setDashStyle = function(el,visProp) {
    var dashStyle = el.visProp['dash'], node = el.rendNode;
    if(el.visProp['dash'] > 0) {
        node.setAttributeNS(null, 'stroke-dasharray', this.dashArray[dashStyle-1]);
    }
    else {
        if(node.hasAttributeNS(null, 'stroke-dasharray')) {
            node.removeAttributeNS(null, 'stroke-dasharray');
        }
    }    
};

JXG.SVGRenderer.prototype.setGridDash = function(id) {
    var node = document.getElementById(id);
    this.setPropertyPrimitive(node,'stroke-dasharray', '5, 5'); 
};

JXG.SVGRenderer.prototype.createPrimitive = function(type,id) {
    var node = this.container.ownerDocument.createElementNS(this.svgNamespace, type);
    node.setAttributeNS(null, 'id', id);
    node.style.position = 'absolute';
    if (type=='path') {
        node.setAttributeNS(null, 'stroke-linecap', 'butt');
        node.setAttributeNS(null, 'stroke-linejoin', 'round');
        //node.setAttributeNS(null, 'shape-rendering', 'geometricPrecision'); // 'crispEdges'
    }
    return node;
};

JXG.SVGRenderer.prototype.createArrowHead = function(el,idAppendix) {
    var id = el.id+'Triangle',
        node2, node3;
        
    if (idAppendix!=null) { id += idAppendix; }
    node2 = this.createPrimitive('marker',id);
    node2.setAttributeNS(null, 'viewBox', '0 0 10 6');
    node2.setAttributeNS(null, 'refY', '3');
    node2.setAttributeNS(null, 'markerUnits', 'strokeWidth');
    node2.setAttributeNS(null, 'markerHeight', '6');
    node2.setAttributeNS(null, 'markerWidth', '6');
    node2.setAttributeNS(null, 'orient', 'auto');
    node2.setAttributeNS(null, 'stroke', el.visProp['strokeColor']);
    node2.setAttributeNS(null, 'stroke-opacity', el.visProp['strokeOpacity']);            
    node2.setAttributeNS(null, 'fill', el.visProp['strokeColor']);
    node2.setAttributeNS(null, 'fill-opacity', el.visProp['strokeOpacity']);    
    node3 = this.container.ownerDocument.createElementNS(this.svgNamespace,'path');
    if (idAppendix=='End') {
        node2.setAttributeNS(null, 'refX', '0');
        node3.setAttributeNS(null, 'd', 'M 0 3 L 10 6 L 10 0 z');
    } else {
        node2.setAttributeNS(null, 'refX', '10');
        node3.setAttributeNS(null, 'd', 'M 0 0 L 10 3 L 0 6 z');
    }
    node2.appendChild(node3);
    return node2;
};

JXG.SVGRenderer.prototype.makeArrow = function(node,el,idAppendix) {
    var node2 = this.createArrowHead(el,idAppendix);
    this.defs.appendChild(node2);
    node.setAttributeNS(null, 'marker-end', 'url(#'+el.id+'Triangle)');
    el.rendNodeTriangle = node2;
};

JXG.SVGRenderer.prototype.makeArrows = function(el) {
    var node2;
    if (el.visPropOld['firstArrow']==el.visProp['firstArrow'] && el.visPropOld['lastArrow']==el.visProp['lastArrow']) {
        return;
    }
    if(el.visProp['firstArrow']) {
        node2 = el.rendNodeTriangleStart;
        if(node2 == null) {
            node2 = this.createArrowHead(el,'End');
            this.defs.appendChild(node2);            
            el.rendNodeTriangleStart = node2;
            el.rendNode.setAttributeNS(null, 'marker-start', 'url(#'+el.id+'TriangleEnd)');    
        }    
    }
    else {
        node2 = el.rendNodeTriangleStart;
        if(node2 != null) {
            this.remove(node2);
        }
    }
    if(el.visProp['lastArrow']) {
        node2 = el.rendNodeTriangleEnd;
        if(node2 == null) {
            node2 = this.createArrowHead(el,'Start');
            this.defs.appendChild(node2);            
            el.rendNodeTriangleEnd = node2;
            el.rendNode.setAttributeNS(null, 'marker-end', 'url(#'+el.id+'TriangleStart)'); 
        }    
    }
    else {
        node2 = el.rendNodeTriangleEnd;
        if(node2 != null) {
            this.remove(node2);
        }        
    }
    el.visPropOld['firstArrow'] = el.visProp['firstArrow'];
    el.visPropOld['lastArrow'] = el.visProp['lastArrow'];
};

JXG.SVGRenderer.prototype.updateLinePrimitive = function(node,p1x,p1y,p2x,p2y) {
    node.setAttributeNS(null, 'x1', p1x);
    node.setAttributeNS(null, 'y1', p1y);
    node.setAttributeNS(null, 'x2', p2x);
    node.setAttributeNS(null, 'y2', p2y);    
};

JXG.SVGRenderer.prototype.updateCirclePrimitive = function(node,x,y,r) {
    node.setAttributeNS(null, 'cx', (x));
    node.setAttributeNS(null, 'cy', (y));
    node.setAttributeNS(null, 'r', (r));
};

JXG.SVGRenderer.prototype.updateEllipsePrimitive = function(node,x,y,rx,ry) {
    node.setAttributeNS(null, 'cx', (x));
    node.setAttributeNS(null, 'cy', (y));
    node.setAttributeNS(null, 'rx', (rx));
    node.setAttributeNS(null, 'ry', (ry));
};

JXG.SVGRenderer.prototype.updateRectPrimitive = function(node,x,y,w,h) {
    node.setAttributeNS(null, 'x', (x));
    node.setAttributeNS(null, 'y', (y));
    node.setAttributeNS(null, 'width', (w));
    node.setAttributeNS(null, 'height', (h));
};

JXG.SVGRenderer.prototype.updatePathPrimitive = function(node, pointString, board) {  // board not necessary in SVG
    /*
    node.setAttributeNS(null, 'stroke-linecap', 'butt');
    node.setAttributeNS(null, 'stroke-linejoin', 'round');
    //node.setAttributeNS(null, 'shape-rendering', 'geometricPrecision');
    //node.setAttributeNS(null, 'shape-rendering', 'crispEdges');
    */
    node.setAttributeNS(null, 'd', pointString);
};

JXG.SVGRenderer.prototype.updatePathStringPrimitive = function(el) {
    var symbm = ' M ',
        symbl = ' L ',
        nextSymb = symbm,
        maxSize = 5000.0,
        pStr = '',
        //h = 3*el.board.canvasHeight,
        //w = 100*el.board.canvasWidth,
        i, scr, 
        isNoPlot = (el.curveType!='plot'),
        //isFunctionGraph = (el.curveType=='functiongraph'),
        len;

    if (el.numberPoints<=0) { return ''; }
    
    if (isNoPlot && el.board.options.curve.RDPsmoothing) {
        el.points = this.RamenDouglasPeuker(el.points,0.5);
    }
    len = Math.min(el.points.length,el.numberPoints);
    for (i=0; i<len; i++) {
        scr = el.points[i].scrCoords;
        //if (isNaN(scr[1]) || isNaN(scr[2]) /*|| Math.abs(scr[1])>w || (isFunctionGraph && (scr[2]>h || scr[2]<-0.5*h))*/ ) {  // PenUp
        if (isNaN(scr[1]) || isNaN(scr[2])) {  // PenUp
            nextSymb = symbm;
        } else {
            // Chrome has problems with values  being too far away.
            if (scr[1]>maxSize) { scr[1] = maxSize; }
            else if (scr[1]<-maxSize) { scr[1] = -maxSize; }
            if (scr[2]>maxSize) { scr[2] = maxSize; }
            else if (scr[2]<-maxSize) { scr[2] = -maxSize; }
            
            pStr += [nextSymb,scr[1],' ',scr[2]].join(''); // Attention: first coordinate may be inaccurate if far way
            nextSymb = symbl;
        }
    }
    return pStr;
};

JXG.SVGRenderer.prototype.updatePathStringPoint = function(el, size, type) {
    var s = '';
    if(type == 'x') {
        s = 'M ' + (el.coords.scrCoords[1]-size) + ' ' + (el.coords.scrCoords[2]-size) + ' L ' + 
        (el.coords.scrCoords[1]+size) + ' ' + (el.coords.scrCoords[2]+size) + ' M ' + 
        (el.coords.scrCoords[1]+size) + ' ' + (el.coords.scrCoords[2]-size) + ' L ' +
        (el.coords.scrCoords[1]-size) + ' ' + (el.coords.scrCoords[2]+size);
    }
    else if(type == '+') {
        s = 'M ' + (el.coords.scrCoords[1]-size) + ' ' + (el.coords.scrCoords[2]) + ' L ' + 
        (el.coords.scrCoords[1]+size) + ' ' + (el.coords.scrCoords[2]) + ' M ' + 
        (el.coords.scrCoords[1]) + ' ' + (el.coords.scrCoords[2]-size) + ' L ' +
        (el.coords.scrCoords[1]) + ' ' + (el.coords.scrCoords[2]+size);    
    }
    else if(type == 'diamond') {
        s = 'M ' + (el.coords.scrCoords[1]-size) + ' ' + (el.coords.scrCoords[2]) + ' L ' + 
        (el.coords.scrCoords[1]) + ' ' + (el.coords.scrCoords[2]+size) + ' L ' + 
        (el.coords.scrCoords[1]+size) + ' ' + (el.coords.scrCoords[2]) + ' L ' +
        (el.coords.scrCoords[1]) + ' ' + (el.coords.scrCoords[2]-size) + ' Z ';
    }
    else if(type == 'A') {
        s = 'M ' + (el.coords.scrCoords[1]) + ' ' + (el.coords.scrCoords[2]-size) + ' L ' + 
        (el.coords.scrCoords[1]-size*Math.sqrt(3)/2) + ' ' + (el.coords.scrCoords[2]+size/2) + ' L ' + 
        (el.coords.scrCoords[1]+size*Math.sqrt(3)/2) + ' ' + (el.coords.scrCoords[2]+size/2) + ' Z ';
    } 
    else if(type == 'v') {
        s = 'M ' + (el.coords.scrCoords[1]) + ' ' + (el.coords.scrCoords[2]+size) + ' L ' + 
        (el.coords.scrCoords[1]-size*Math.sqrt(3)/2) + ' ' + (el.coords.scrCoords[2]-size/2) + ' L ' + 
        (el.coords.scrCoords[1]+size*Math.sqrt(3)/2) + ' ' + (el.coords.scrCoords[2]-size/2) + ' Z ';
    }   
    else if(type == '>') {
        s = 'M ' + (el.coords.scrCoords[1]+size) + ' ' + (el.coords.scrCoords[2]) + ' L ' + 
        (el.coords.scrCoords[1]-size/2) + ' ' + (el.coords.scrCoords[2]-size*Math.sqrt(3)/2) + ' L ' + 
        (el.coords.scrCoords[1]-size/2) + ' ' + (el.coords.scrCoords[2]+size*Math.sqrt(3)/2) + ' Z ';
    }
    else if(type == '<') {
        s = 'M ' + (el.coords.scrCoords[1]-size) + ' ' + (el.coords.scrCoords[2]) + ' L ' + 
        (el.coords.scrCoords[1]+size/2) + ' ' + (el.coords.scrCoords[2]-size*Math.sqrt(3)/2) + ' L ' + 
        (el.coords.scrCoords[1]+size/2) + ' ' + (el.coords.scrCoords[2]+size*Math.sqrt(3)/2) + ' Z ';
    }
    return s;
}

JXG.SVGRenderer.prototype.updatePolygonePrimitive = function(node, el) {
    var pStr = '', 
        scrCoords, i,
        len = el.vertices.length;
        
    node.setAttributeNS(null, 'stroke', 'none');
    for(i=0; i<len-1; i++) {
        scrCoords = el.vertices[i].coords.scrCoords;
        pStr = pStr + scrCoords[1] + "," + scrCoords[2];
        if(i<len-2) { pStr += " "; }
    }
    node.setAttributeNS(null, 'points', pStr);
};

JXG.SVGRenderer.prototype.appendChildPrimitive = function(node,level) {
    if (typeof level=='undefined') { // trace nodes have level not set
        level = 0;                         
    } else if (level>=JXG.Options.layer.numlayers) { 
        level = JXG.Options.layer.numlayers-1;
    }
    this.layer[level].appendChild(node);
};

JXG.SVGRenderer.prototype.setPropertyPrimitive = function(node,key,val) {
    if (key=='stroked') {
        return;
    }
    node.setAttributeNS(null, key, val);
};

JXG.SVGRenderer.prototype.drawVerticalGrid = function(topLeft, bottomRight, gx, board) {
    var node = this.createPrimitive('path', 'gridx'),
        gridArr = '';
        
    while(topLeft.scrCoords[1] < bottomRight.scrCoords[1] + gx - 1) { 
        gridArr += ' M ' + topLeft.scrCoords[1] + ' ' + 0 + ' L ' + topLeft.scrCoords[1] + ' ' + board.canvasHeight+' ';
        topLeft.setCoordinates(JXG.COORDS_BY_SCREEN, [topLeft.scrCoords[1] + gx, topLeft.scrCoords[2]]);   
    }
    this.updatePathPrimitive(node, gridArr, board);
    return node;
};

JXG.SVGRenderer.prototype.drawHorizontalGrid = function(topLeft, bottomRight, gy, board) {
    var node = this.createPrimitive('path', 'gridy'),
        gridArr = '';
        
    while(topLeft.scrCoords[2] <= bottomRight.scrCoords[2] + gy - 1) {
        gridArr += ' M ' + 0 + ' ' + topLeft.scrCoords[2] + ' L ' + board.canvasWidth + ' ' + topLeft.scrCoords[2]+' ';
        topLeft.setCoordinates(JXG.COORDS_BY_SCREEN, [topLeft.scrCoords[1], topLeft.scrCoords[2] + gy]);
    }
    this.updatePathPrimitive(node, gridArr, board);
    return node;
};

JXG.SVGRenderer.prototype.appendNodesToElement = function(element, type) {
    element.rendNode = document.getElementById(element.id);
};


/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
--------------------------------------------------------------------
von AbstractRenderer abgeleitete Zeichenklasse
fuer Browser mit VML-Elementen (Internet Explorer)
--------------------------------------------------------------------
*/

JXG.VMLRenderer = function(container) {
    this.constructor();
    
    this.container = container;
    this.container.style.overflow = 'hidden';
    this.container.onselectstart = function () { return false; };
    
    this.resolution = 10; // Paths are drawn with a a resolution of this.resolution/pixel.
  
    // Add VML includes and namespace
    // Original: IE <=7
    container.ownerDocument.namespaces.add("jxgvml", "urn:schemas-microsoft-com:vml");
    //container.ownerDocument.createStyleSheet().addRule("v\\:*", "behavior: url(#default#VML);");

    this.container.ownerDocument.createStyleSheet().addRule(".jxgvml", "behavior:url(#default#VML)");
    try {
        !container.ownerDocument.namespaces.jxgvml && container.ownerDocument.namespaces.add("jxgvml", "urn:schemas-microsoft-com:vml");
        this.createNode = function (tagName) {
            return container.ownerDocument.createElement('<jxgvml:' + tagName + ' class="jxgvml">');
        };
    } catch (e) {
        this.createNode = function (tagName) {
            return container.ownerDocument.createElement('<' + tagName + ' xmlns="urn:schemas-microsoft.com:vml" class="jxgvml">');
        };
    }
        
    // um Dashes zu realisieren
    this.dashArray = ['Solid', '1 1', 'ShortDash', 'Dash', 'LongDash', 'ShortDashDot', 'LongDashDot'];    
};

JXG.VMLRenderer.prototype = new JXG.AbstractRenderer;

JXG.VMLRenderer.prototype.setAttr = function(node, key, val, val2) {
    try {
        if (document.documentMode==8) {
            node[key] = val;
        } else {
            node.setAttribute(key,val,val2);
        }
    } catch (e) {
        //document.getElementById('debug').innerHTML += node.id+' '+key+' '+val+'<br>\n';
    }
};

JXG.VMLRenderer.prototype.setShadow = function(el) {
    var nodeShadow = el.rendNodeShadow;
    
    if (!nodeShadow) return;                          // Added 29.9.09. A.W.
    if (el.visPropOld['shadow']==el.visProp['shadow']) {
        return;
    }
    if(el.visProp['shadow']) {
        this.setAttr(nodeShadow, 'On', 'True');
        this.setAttr(nodeShadow, 'Offset', '3pt,3pt');
        this.setAttr(nodeShadow, 'Opacity', '60%');
        this.setAttr(nodeShadow, 'Color', '#aaaaaa');
    }
    else {
        this.setAttr(nodeShadow, 'On', 'False');
    }
    el.visPropOld['shadow']=el.visProp['shadow'];
};

JXG.VMLRenderer.prototype.setGradient = function(el) {
    var nodeFill = el.rendNodeFill;
    if(el.type == JXG.OBJECT_TYPE_ARC || el.type == JXG.OBJECT_TYPE_ANGLE) {
        nodeFill = el.rendNode2Fill;
    }    
    
    if(el.visProp['gradient'] == 'linear') {
        this.setAttr(nodeFill, 'type', 'gradient');
        this.setAttr(nodeFill, 'color2', el.visProp['gradientSecondColor']);
        this.setAttr(nodeFill, 'opacity2', el.visProp['gradientSecondOpacity']);
        this.setAttr(nodeFill, 'angle', el.visProp['gradientAngle']);
    }
    else if (el.visProp['gradient'] == 'radial') {
        this.setAttr(nodeFill, 'type','gradientradial');
        this.setAttr(nodeFill, 'color2',el.visProp['gradientSecondColor']);
        this.setAttr(nodeFill, 'opacity2',el.visProp['gradientSecondOpacity']);
        this.setAttr(nodeFill, 'focusposition', el.visProp['gradientPositionX']*100+'%,'+el.visProp['gradientPositionY']*100+'%');
        this.setAttr(nodeFill, 'focussize', '0,0');
    }
    else {
        this.setAttr(nodeFill, 'type','solid');
    }
};

JXG.VMLRenderer.prototype.updateGradient = function(el) {}; // Not needed in VML;

JXG.VMLRenderer.prototype.addShadowToGroup = function(groupname, board) {
    var el, pEl;
    if(groupname == "lines") {
        for(el in board.objects) {
            pEl = board.objects[el];
            if(pEl.elementClass == JXG.OBJECT_CLASS_LINE) {
                this.addShadowToElement(pEl);
            }
        }
    }
    else if(groupname == "points") {
        for(el in board.objects) {
            pEl = board.objects[el];
            if(pEl.elementClass == JXG.OBJECT_CLASS_POINT) {
                this.addShadowToElement(pEl);
            }
        }
    }
    else if(groupname == "circles") {
        for(el in board.objects) {
            pEl = board.objects[el];
            if(pEl.elementClass == JXG.OBJECT_CLASS_CIRCLE) {
                this.addShadowToElement(pEl);
            }
        }
    }    
    board.fullUpdate();
};

JXG.VMLRenderer.prototype.displayCopyright = function(str,fontsize) {
    var node, t;
    
    //node = this.container.ownerDocument.createElement('v:textbox');
    node = this.createNode('textbox');
    node.style.position = 'absolute';
    this.setAttr(node,'id', 'licenseText');
    
    node.style.left = 20;
    node.style.top = (2);
    node.style.fontSize = (fontsize);
    node.style.color = '#356AA0';
    node.style.fontFamily = 'Arial,Helvetica,sans-serif';
    this.setAttr(node,'opacity','30%');
    node.style.filter = 'alpha(opacity = 30)';
    
    t = document.createTextNode(str);
    node.appendChild(t);
    this.appendChildPrimitive(node,0);
};
JXG.VMLRenderer.prototype.drawInternalText = function(el) {
    var node;
    node = this.createNode('textbox');
    node.style.position = 'absolute';
    if (document.documentMode==8) {    
        node.setAttribute('class', 'JXGtext');
    } else {
        node.setAttribute('className', 9);
    }
    el.rendNodeText = document.createTextNode('');
    node.appendChild(el.rendNodeText);
    this.appendChildPrimitive(node,9);
    return node;
};

JXG.VMLRenderer.prototype.updateInternalText = function(/** JXG.Text */ el) { 
    el.rendNode.style.left = (el.coords.scrCoords[1])+'px'; 
    el.rendNode.style.top = (el.coords.scrCoords[2] - this.vOffsetText)+'px'; 
    el.updateText();
    if (el.htmlStr!= el.plaintextStr) {
        el.rendNodeText.data = el.plaintextStr;
        el.htmlStr = el.plaintextStr;
    }
};

JXG.VMLRenderer.prototype.drawTicks = function(ticks) {
    var ticksNode = this.createPrimitive('path', ticks.id);
    this.appendChildPrimitive(ticksNode,ticks.layer);
    //ticks.rendNode = ticksNode;
    this.appendNodesToElement(ticks, 'path');
};

JXG.VMLRenderer.prototype.updateTicks = function(axis,dxMaj,dyMaj,dxMin,dyMin) {
    var tickArr = [], i, len, c, ticks;
    
    len = axis.ticks.length;
    for (i=0; i<len; i++) {
        c = axis.ticks[i];
        if(c.major) {
            if (axis.labels[i].visProp['visible']) this.drawText(axis.labels[i]);        
            tickArr.push(' m ' + Math.round(this.resolution*(c.scrCoords[1]+dxMaj)) + 
                         ', ' + Math.round(this.resolution*(c.scrCoords[2]-dyMaj)) + 
                         ' l ' + Math.round(this.resolution*(c.scrCoords[1]-dxMaj)) + 
                         ', ' + Math.round(this.resolution*(c.scrCoords[2]+dyMaj))+' ');
        }
        else
            tickArr.push(' m ' + Math.round(this.resolution*(c.scrCoords[1]+dxMin)) + 
                         ', ' + Math.round(this.resolution*(c.scrCoords[2]-dyMin)) + 
                         ' l ' + Math.round(this.resolution*(c.scrCoords[1]-dxMin)) + 
                         ', ' + Math.round(this.resolution*(c.scrCoords[2]+dyMin))+' ');
    }

    ticks = document.getElementById(axis.id);
    if(ticks == null) {
        ticks = this.createPrimitive('path', axis.id);
        this.appendChildPrimitive(ticks,axis.layer);
        this.appendNodesToElement(axis,'path');
    } 
    this.setAttr(ticks,'stroked', 'true');
    this.setAttr(ticks,'strokecolor', axis.visProp['strokeColor'], 1);
    this.setAttr(ticks,'strokeweight', axis.visProp['strokeWidth']);   
    //ticks.setAttributeNS(null, 'stroke-opacity', axis.visProp['strokeOpacity']);
    this.updatePathPrimitive(ticks, tickArr, axis.board);
};

JXG.VMLRenderer.prototype.drawArcLine = function(id, radius, angle1, angle2, midpoint, board, el) {
    //var node = this.createPrimitive('arc',id);  doesn't work here
    var node = this.createNode('arc'),
        fillNode = this.createNode('fill'),
        strokeNode = this.createNode('stroke'),
        shadowNode = this.createNode('shadow');

    this.setAttr(node,'id', id);
    this.setAttr(fillNode,'id', id+'_fill');
    this.setAttr(strokeNode,'id', id+'_stroke');
    this.setAttr(shadowNode,'id', id+'_shadow');

    node.appendChild(fillNode);
    node.appendChild(strokeNode);
    node.appendChild(shadowNode);   
    
    el.rendNode = node;
    el.rendNodeFill = fillNode;
    el.rendNodeStroke = strokeNode;
    el.rendNodeShadow = shadowNode;    
    
    node.style.position = 'absolute';
    this.setAttr(node,'filled', 'false');

    node.style.left = (midpoint.coords.scrCoords[1] - Math.round(radius * board.stretchX)) + 'px'; 
    node.style.top = (midpoint.coords.scrCoords[2] - Math.round(radius * board.stretchY))  + 'px'; 
    node.style.width = (Math.round(radius * board.stretchX)*2) + 'px'; 
    node.style.height = (Math.round(radius * board.stretchY)*2) + 'px';  
    this.setAttr(node,'startangle', angle1);
    this.setAttr(node,'endangle', angle2);  
    
    return node;
};

JXG.VMLRenderer.prototype.drawArcFill = function(id, radius, midpoint, point2, point3, board, element) {
    // createPrimitive doesn't work here...
    var x, y, pathString,
        pathNode = this.createNode('path'),
        node2 = this.createNode('shape'),
        fillNode = this.createNode('fill'),
        strokeNode = this.createNode('stroke'),
        shadowNode = this.createNode('shadow');
 
    id = id+'sector';
    this.setAttr(pathNode, 'id', id+'_path');        
    this.setAttr(fillNode,'id', id+'_fill');
    this.setAttr(strokeNode,'id', id+'_stroke');
    this.setAttr(shadowNode,'id', id+'_shadow');    
    this.setAttr(node2,'id', id);
 
    node2.appendChild(fillNode);
    node2.appendChild(strokeNode);
    node2.appendChild(shadowNode);   
    node2.appendChild(pathNode);
    node2.style.position = 'absolute';    
    
    element.rendNode2 = node2;
    element.rendNode2Fill = fillNode;
    element.rendNode2Stroke = strokeNode;
    element.rendNode2Shadow = shadowNode;
    element.rendNode2Path = pathNode;
        
    this.setAttr(node2,'stroked', 'false');

    x = Math.round(radius * board.stretchX); // Breite des umgebenden Rechtecks?
    y = Math.round(radius * board.stretchY); // Hoehe des umgebenden Rechtecks?

    node2.style.width = x;
    node2.style.height = y;
    this.setAttr(node2,'coordsize', x+','+y);

    pathString = 'm ' + midpoint.coords.scrCoords[1] + ',' + midpoint.coords.scrCoords[2] + ' l ';  
    pathString += point2.coords.scrCoords[1] + ',' + point2.coords.scrCoords[2] + ' at ';
    pathString += (midpoint.coords.scrCoords[1]-x) + ',' + (midpoint.coords.scrCoords[2]-y) + ',';
    pathString += (midpoint.coords.scrCoords[1]+x) + ',' + (midpoint.coords.scrCoords[2]+y);
    pathString += ' ' + point2.coords.scrCoords[1] + ',' + point2.coords.scrCoords[2];
    pathString += ', ' + point3.coords.scrCoords[1] + ',' + point3.coords.scrCoords[2] + ' l ';
    pathString += midpoint.coords.scrCoords[1] + ',' + midpoint.coords.scrCoords[2] + ' x e';
    
    this.setAttr(pathNode,'v', pathString);
    
    return node2;
};

JXG.VMLRenderer.prototype.drawArc = function(el) { 
    var radius, p = {}, angle1, angle2, node, nodeStroke, node2, p4 = {};

    JXG.clearVisPropOld(el);

    /* some computations */
    radius = el.Radius();  
    p.coords = new JXG.Coords(JXG.COORDS_BY_USER, 
                          [el.midpoint.coords.usrCoords[1], el.board.origin.scrCoords[2]/el.board.stretchY],
                          el.board);
    angle2 = el.board.algebra.trueAngle(el.point2, el.midpoint, p);
    angle1 = el.board.algebra.trueAngle(el.point3, el.midpoint, p);
    if(angle2 < angle1) {
        angle1 -= 360;
    }
    
    /* arc line */
    node = this.drawArcLine(el.id, radius, angle1, angle2, el.midpoint, el.board, el);
    
    /* arrows at the ends of the arc line */
    nodeStroke = el.rendNodeStroke;
    if(el.visProp['lastArrow']) {        
        this.setAttr(nodeStroke,'endarrow', 'block');
        this.setAttr(nodeStroke,'endarrowlength', 'long');
    }
    if(el.visProp['firstArrow']) {        
        this.setAttr(nodeStroke,'startarrow', 'block');
        this.setAttr(nodeStroke,'startarrowlength', 'long');
    }
    
    /* stroke color and width */
    this.setObjectStrokeColor(el,el.visProp['strokeColor'],el.visProp['strokeOpacity']);
    this.setObjectStrokeWidth(el,el.visProp['strokeWidth']);
    
    /* dashstyle and shadow */
    this.setDashStyle(el,el.visProp);  
    this.setShadow(el); 
   
    /* arc fill */
    p4.coords = el.board.algebra.projectPointToCircle(el.point3,el);      
    node2 = this.drawArcFill(el.id, radius, el.midpoint, el.point2, p4, el.board, el);
    
    /* fill props */
    this.setObjectFillColor(el, el.visProp['fillColor'], el.visProp['fillOpacity']);
    this.setGradient(el);    
    
    /* append nodes */
    this.appendChildPrimitive(node,el.layer); //arc
    this.appendChildPrimitive(node2,el.layer); //fill
    
    /* draft mode */
    if(el.visProp['draft']) {
       this.setDraft(el);
    }
    if(!el.visProp['visible']) {
        el.hideElement(el);
    }
};

/**
 * Updates properties of an arc that already exists.
 * @param {JXG.Arc} arc Reference to an arc object, that has to be updated.
 * @see JXG.Arc
 * @see #drawArc
 */
JXG.AbstractRenderer.prototype.updateArc = function(el) { 
    // AW: brutaler fix der update-Methode...
    this.remove(el.rendNode);     
    this.remove(el.rendNode2);
    this.drawArc(el);
    return;
};

JXG.VMLRenderer.prototype.drawAngle = function(el) {
    var circle  = {}, projectedP1, projectedP3, p = {}, 
        angle1, angle2, node, tmp, nodeStroke,
        p1 = {}, p3 = {}, node2;

    JXG.clearVisPropOld(el);
    
    /* some computations */
    // um projectToCircle benutzen zu koennen...
    circle.midpoint = el.point2;
    circle.Radius = function() {
        return el.radius;
    };
    //-----------------
    // deprecated:
    circle.getRadius = function() {
        return el.radius;
    };
    //-----------------
    projectedP1 = el.board.algebra.projectPointToCircle(el.point1,circle);
    projectedP3 = el.board.algebra.projectPointToCircle(el.point3,circle);  
    
    p.coords = new JXG.Coords(JXG.COORDS_BY_USER, 
                          [el.point2.coords.usrCoords[1], el.board.origin.scrCoords[2]/(el.board.stretchY)],
                          el.board);
    angle2 = el.board.algebra.trueAngle(el.point1, el.point2, p);
    angle1 = el.board.algebra.trueAngle(el.point3, el.point2, p);
    if(angle2 < angle1) {
        angle1 -= 360;
    }    

    /* arc line */
    node = this.drawArcLine(el.id, el.radius, angle1, angle2, el.point2, el.board, el);

    /* stroke color and width */
    this.setObjectStrokeColor(el,el.visProp['strokeColor'],el.visProp['strokeOpacity']);
    this.setObjectStrokeWidth(el,el.visProp['strokeWidth']);
    
    /* dashstyle and shadow */
    tmp = el.visProp['dash'];
    nodeStroke = el.rendNodeStroke;    
    this.setAttr(nodeStroke,'dashstyle', this.dashArray[tmp]);     
    this.setShadow(el); 
   
    /* arc fill */
    p1.coords = projectedP1;  
    p3.coords = projectedP3;
    node2 = this.drawArcFill(el.id, el.radius, el.point2, p1, p3, el.board, el);   

    /* fill props */
    this.setObjectFillColor(el, el.visProp['fillColor'], el.visProp['fillOpacity']);
    
    /* append nodes */
    this.appendChildPrimitive(node,el.layer); //arc
    this.appendChildPrimitive(node2,el.layer); //fill
    
    /* draft mode */
    if(el.visProp['draft']) {
       this.setDraft(el);
    }

    if(!el.visProp['visible']) {
        el.hideElement(el);
    }
};

JXG.VMLRenderer.prototype.updateAngle = function(el) {
    // erstmal nur der brutale Weg... 
    this.remove(el.rendNode);
    this.remove(el.rendNode2);  
    this.drawAngle(el);
    return;
};

JXG.VMLRenderer.prototype.drawImage = function(el) {
    // IE 8: Bilder ueber data URIs werden bis 32kB unterstuetzt.
    var node, url = el.url; //'data:image/png;base64,' + el.imageBase64String;    
    
    node = this.container.ownerDocument.createElement('img');
    node.style.position = 'absolute';
    this.setAttr(node,'id', el.id);

    this.setAttr(node,'src',url);
    this.container.appendChild(node);
    this.appendChildPrimitive(node,el.layer);
    node.style.filter = "progid:DXImageTransform.Microsoft.Matrix(M11='1.0', sizingMethod='auto expand')";
    el.rendNode = node;
    this.updateImage(el);
};

JXG.VMLRenderer.prototype.transformImage = function(el,t) {
    var node = el.rendNode, 
        m;
    m = this.joinTransforms(el,t);
    node.style.left = (el.coords.scrCoords[1] + m[1][0]) + 'px'; 
    node.style.top = (el.coords.scrCoords[2]-el.size[1] + m[2][0]) + 'px';    
    node.filters.item(0).M11 = m[1][1];
    node.filters.item(0).M12 = m[1][2];
    node.filters.item(0).M21 = m[2][1];
    node.filters.item(0).M22 = m[2][2];
};

JXG.VMLRenderer.prototype.joinTransforms = function(el,t) {
    var m = [[1,0,0],[0,1,0],[0,0,1]], 
        i,
        len = t.length;
        
    for (i=0;i<len;i++) {
        m = JXG.Math.matMatMult(t[i].matrix,m);
    }
    return m;
};

JXG.VMLRenderer.prototype.transformImageParent = function(el,m) {};

/*
JXG.VMLRenderer.prototype.removeGrid = function(board) { 
    var c = document.getElementById('gridx');
    this.remove(c);

    c = document.getElementById('gridy');
    this.remove(c);

    board.hasGrid = false;
};
*/

JXG.VMLRenderer.prototype.hide = function(el) {
    var node = el.rendNode;
    node.style.visibility = "hidden"; 
    if(el.type == JXG.OBJECT_TYPE_ARC || el.type == JXG.OBJECT_TYPE_ANGLE) {
        node = el.rendNode2; 
        node.style.visibility = "hidden";         
    }
};

JXG.VMLRenderer.prototype.show = function(el) {
    var node = el.rendNode;
    node.style.visibility = "inherit";  
    if(el.type == JXG.OBJECT_TYPE_ARC || el.type == JXG.OBJECT_TYPE_ANGLE) {
        node = el.rendNode2; 
        node.style.visibility = "inherit";         
    }
};

JXG.VMLRenderer.prototype.setDashStyle = function(el,visProp) {
    var node;
    if(visProp['dash'] >= 0) {
        node = el.rendNodeStroke;
        this.setAttr(node,'dashstyle', this.dashArray[visProp['dash']]);
    }
};
 
JXG.VMLRenderer.prototype.setObjectStrokeColor = function(el, color, opacity) {
    var c = this.eval(color), 
        o = this.eval(opacity), 
        node, nodeStroke;

    o = (o>0)?o:0;

    if (el.visPropOld['strokeColor']==c && el.visPropOld['strokeOpacity']==o) {
        return;
    }
    if(el.type == JXG.OBJECT_TYPE_TEXT) {
        el.rendNode.style.color = c;
    }        
    else {       
        node = el.rendNode;
        this.setAttr(node,'stroked', 'true');
        this.setAttr(node,'strokecolor', c);
        
        if(el.id == 'gridx') {
            nodeStroke = document.getElementById('gridx_stroke')
        }
        else if(el.id == 'gridy') {
            nodeStroke = document.getElementById('gridy_stroke')
        }
        else {
            nodeStroke = el.rendNodeStroke;
        }
        if (o!=undefined) {
            this.setAttr(nodeStroke,'opacity', (o*100)+'%');  
            
        }
    }
    el.visPropOld['strokeColor'] = c;
    el.visPropOld['strokeOpacity'] = o;
};

JXG.VMLRenderer.prototype.setObjectFillColor = function(el, color, opacity) {
    var c = this.eval(color), 
        o = this.eval(opacity);

    o = (o>0)?o:0;

    if (el.visPropOld['fillColor']==c && el.visPropOld['fillOpacity']==o) {
        return;
    }
    
    if(el.type == JXG.OBJECT_TYPE_ARC || el.type == JXG.OBJECT_TYPE_ANGLE) {
        if(c == 'none') {
             this.setAttr(el.rendNode2,'filled', 'false');
        }
        else {
            this.setAttr(el.rendNode2,'filled', 'true');
            this.setAttr(el.rendNode2,'fillcolor', c); 
            if (o!=undefined) {
                 this.setAttr(el.rendNode2Fill,'opacity', (o*100)+'%');
            }
        }
    }
    else {
        if(c == 'none') {
            this.setAttr(el.rendNode,'filled', 'false');
        }
        else {
            this.setAttr(el.rendNode,'filled', 'true');
            this.setAttr(el.rendNode,'fillcolor', c); 
            if (o!=undefined && el.rendNodeFill) {  // Added el.rendNodeFill 29.9.09  A.W.
                this.setAttr(el.rendNodeFill,'opacity', (o*100)+'%');
            }
        }
    }
    el.visPropOld['fillColor'] = c;
    el.visPropOld['fillOpacity'] = o;
};

JXG.VMLRenderer.prototype.remove = function(node) {
  if (node!=null) node.removeNode(true);
};

JXG.VMLRenderer.prototype.suspendRedraw = function() {
    this.container.style.display='none';
};

JXG.VMLRenderer.prototype.unsuspendRedraw = function() {
    this.container.style.display='';
};

JXG.VMLRenderer.prototype.setAttributes = function(node,props,vmlprops,visProp) {
    var val, i, p
        len = props.length;

    for (i=0;i<len;i++) {
        p = props[i];
        if (visProp[p]!=null) {
            val = this.eval(visProp[p]);
            val = (val>0)?val:0;
            this.setAttr(node,vmlprops[i], val);
        }
    }
};

JXG.VMLRenderer.prototype.setGridDash = function(id, node) {
    var node = document.getElementById(id+'_stroke');
    this.setAttr(node,'dashstyle', 'Dash');
};

/**
 * Sets an elements stroke width.
 * @param {Object} el Reference to the geometry element.
 * @param {int} width The new stroke width to be assigned to the element.
 */
JXG.VMLRenderer.prototype.setObjectStrokeWidth = function(el, width) {
    var w = this.eval(width), 
        node;
    //w = (w>0)?w:0;
    
    if (el.visPropOld['strokeWidth']==w) {
        return;
    }
    
    node = el.rendNode;
    this.setPropertyPrimitive(node,'stroked', 'true');
    if (w!=null) { 
        this.setPropertyPrimitive(node,'stroke-width',w); 
    }
    el.visPropOld['strokeWidth'] = w;
};

JXG.VMLRenderer.prototype.createPrimitive = function(type, id) {
    var node, 
        fillNode = this.createNode('fill'), 
        strokeNode = this.createNode('stroke'), 
        shadowNode = this.createNode('shadow'), 
        pathNode;
    
    this.setAttr(fillNode, 'id', id+'_fill');
    this.setAttr(strokeNode, 'id', id+'_stroke');
    this.setAttr(shadowNode, 'id', id+'_shadow');
    
    if (type=='circle' || type=='ellipse' ) {
        node = this.createNode('oval');
        node.appendChild(fillNode);
        node.appendChild(strokeNode);
        node.appendChild(shadowNode);
    } else if (type == 'polygon' || type == 'path' || type == 'shape' || type == 'line') {    
        node = this.createNode('shape');
        node.appendChild(fillNode);
        node.appendChild(strokeNode);
        node.appendChild(shadowNode);   
        pathNode = this.createNode('path');
        this.setAttr(pathNode, 'id', id+'_path');        
        node.appendChild(pathNode);
    } else {
        node = this.createNode(type);
        node.appendChild(fillNode);
        node.appendChild(strokeNode);
        node.appendChild(shadowNode);
    }
    node.style.position = 'absolute';
    this.setAttr(node, 'id', id);
    
    return node;
};

JXG.VMLRenderer.prototype.appendNodesToElement = function(element, type) {
    if(type == 'shape' || type == 'path' || type == 'polygon') {
        element.rendNodePath = document.getElementById(element.id+'_path');
    }
    element.rendNodeFill = document.getElementById(element.id+'_fill');
    element.rendNodeStroke = document.getElementById(element.id+'_stroke');
    element.rendNodeShadow = document.getElementById(element.id+'_shadow');
    element.rendNode = document.getElementById(element.id);
};

JXG.VMLRenderer.prototype.makeArrow = function(node,el,idAppendix) {
    var nodeStroke = el.rendNodeStroke;
    this.setAttr(nodeStroke, 'endarrow', 'block');
    this.setAttr(nodeStroke, 'endarrowlength', 'long');
};

JXG.VMLRenderer.prototype.makeArrows = function(el) {
    var nodeStroke;
    
    if (el.visPropOld['firstArrow']==el.visProp['firstArrow'] && el.visPropOld['lastArrow']==el.visProp['lastArrow']) {
        return;
    }

    if(el.visProp['firstArrow']) {
        nodeStroke = el.rendNodeStroke;
        this.setAttr(nodeStroke, 'startarrow', 'block');
        this.setAttr(nodeStroke, 'startarrowlength', 'long');                 
    }
    else {
        nodeStroke = el.rendNodeStroke;
        if(nodeStroke != null) {
            this.setAttr(nodeStroke, 'startarrow', 'none');
        }            
    }
    if(el.visProp['lastArrow']) {
        nodeStroke = el.rendNodeStroke;
        this.setAttr(nodeStroke, 'id', el.id+"stroke");
        this.setAttr(nodeStroke, 'endarrow', 'block');
        this.setAttr(nodeStroke, 'endarrowlength', 'long');            
    }
    else {
        nodeStroke = el.rendNodeStroke;
        if(nodeStroke != null) {
            this.setAttr(nodeStroke, 'endarrow', 'none');
        }        
    }    
    el.visPropOld['firstArrow'] = el.visProp['firstArrow'];
    el.visPropOld['lastArrow'] = el.visProp['lastArrow'];
};

JXG.VMLRenderer.prototype.updateLinePrimitive = function(node,p1x,p1y,p2x,p2y,board) {
    /* 
    this.setAttr(node, 'from', [p1x,p1y].join(',')); 
    this.setAttr(node, 'to', [p2x,p2y].join(','));      
    */
    var s, r = this.resolution;
    s = ['m ',r*p1x,', ',r*p1y,' l ',r*p2x,', ',r*p2y];
    this.updatePathPrimitive(node,s,board);
};

JXG.VMLRenderer.prototype.updateCirclePrimitive = function(node,x,y,r) {
    //node.setAttribute('style','left:'+(x-r)+'px; top:'+(y-r)+'px; width:'+(r*2)+'px; height:'+ (r*2)+'px'); 
    node.style.left = (x-r)+'px';
    node.style.top = (y-r)+'px';    
    node.style.width = (r*2)+'px'; 
    node.style.height = (r*2)+'px';   
};

JXG.VMLRenderer.prototype.updateRectPrimitive = function(node,x,y,w,h) {
    node.style.left = (x)+'px';
    node.style.top = (y)+'px';    
    node.style.width = (w)+'px'; 
    node.style.height = (h)+'px';   
};

JXG.VMLRenderer.prototype.updateEllipsePrimitive = function(node,x,y,rx,ry) {
    node.style.left = (x-rx)+'px';
    node.style.top =  (y-ry)+'px'; 
    node.style.width = (rx*2)+'px'; 
    node.style.height = (ry*2)+'px';
};

JXG.VMLRenderer.prototype.updatePathPrimitive = function(node,pointString,board) {
    var x = board.canvasWidth, 
        y = board.canvasHeight;
    node.style.width = x;
    node.style.height = y;
    this.setAttr(node, 'coordsize', [(this.resolution*x),(this.resolution*y)].join(','));
    this.setAttr(node, 'path',pointString.join(""));
};

JXG.VMLRenderer.prototype.updatePathStringPrimitive = function(el) {
    var pStr = [], 
        //h = 3*el.board.canvasHeight, 
        //w = 100*el.board.canvasWidth, 
        i, scr,
        r = this.resolution,
        mround = Math.round,
        symbm = ' m ', 
        symbl = ' l ',
        nextSymb = symbm, 
        isNoPlot = (el.curveType!='plot'),
        //isFunctionGraph = (el.curveType=='functiongraph'),
        len = Math.min(el.numberPoints,8192); // otherwise IE 7 crashes in hilbert.html
    
    if (el.numberPoints<=0) { return ''; }
    if (isNoPlot && el.board.options.curve.RDPsmoothing) {
        el.points = this.RamenDouglasPeuker(el.points,1.0);
    }
    len = Math.min(len,el.points.length);

    for (i=0; i<len; i++) {
        scr = el.points[i].scrCoords;
        if (isNaN(scr[1]) || isNaN(scr[2]) /* || Math.abs(scr[1])>w || (isFunctionGraph && (scr[2]>h || scr[2]<-0.5*h))*/ ) {  // PenUp
            nextSymb = symbm;
        } else {
            // IE has problems with values  being too far away.
            if (scr[1]>20000.0) { scr[1] = 20000.0; }
            else if (scr[1]<-20000.0) { scr[1] = -20000.0; }
            if (scr[2]>20000.0) { scr[2] = 20000.0; }
            else if (scr[2]<-20000.0) { scr[2] = -20000.0; }

            pStr.push([nextSymb,mround(r*scr[1]),', ',mround(r*scr[2])].join(''));
            nextSymb = symbl;
        }
    }
    pStr.push(' e');
    return pStr;
};

JXG.VMLRenderer.prototype.updatePathStringPoint = function(el, size, type) {
    var s = [],
        scr = el.coords.scrCoords,
        r = this.resolution;

    if(type == 'x') {
        s.push(['m ',(r*(scr[1]-size)),', ',(r*(scr[2]-size)),' l ',
        (r*(scr[1]+size)),', ',(r*(scr[2]+size)),' m ',
        (r*(scr[1]+size)),', ',(r*(scr[2]-size)),' l ',
        (r*(scr[1]-size)),', ',(r*(scr[2]+size))].join(''));
    }
    else if(type == '+') {
        s.push(['m ',(r*(scr[1]-size)),', ',(r*(scr[2])),' l ',
        (r*(scr[1]+size)),', ',(r*(scr[2])),' m ',
        (r*(scr[1])),', ',(r*(scr[2]-size)),' l ',
        (r*(scr[1])),', ',(r*(scr[2]+size))].join(''));    
    }
    else if(type == 'diamond') {
        s.push(['m ',(r*(scr[1]-size)),', ',(r*(scr[2])),' l ',
        (r*(scr[1])),', ',(r*(scr[2]+size)),' l ',
        (r*(scr[1]+size)),', ',(r*(scr[2])),' l ',
        (r*(scr[1])),', ',(r*(scr[2]-size)),' x e '
        ].join(''));   
    }
    else if(type == 'A') {
        s.push(['m ',(r*(scr[1])),', ',(r*(scr[2]-size)),' l ',
        Math.round(r*(scr[1]-size*Math.sqrt(3)/2)),', ',(r*(scr[2]+size/2)),' l ',
        Math.round(r*(scr[1]+size*Math.sqrt(3)/2)),', ',(r*(scr[2]+size/2)),' x e '
        ].join(''));           
    } 
    else if(type == 'v') {
        s.push(['m ',(r*(scr[1])),', ',(r*(scr[2]+size)),' l ',
        Math.round(r*(scr[1]-size*Math.sqrt(3)/2)),', ',(r*(scr[2]-size/2)),' l ',
        Math.round(r*(scr[1]+size*Math.sqrt(3)/2)),', ',(r*(scr[2]-size/2)),' x e '
        ].join(''));       
    }   
    else if(type == '>') {
        s.push(['m ',(r*(scr[1]+size)),', ',(r*(scr[2])),' l ',
        (r*(scr[1]-size/2)),', ',Math.round(r*(scr[2]-size*Math.sqrt(3)/2)),' l ',
        (r*(scr[1]-size/2)),', ',Math.round(r*(scr[2]+size*Math.sqrt(3)/2)),
        //' x e '
        ' l ',(r*(scr[1]+size)),', ',(r*(scr[2])) 
        ].join(''));        
    }
    else if(type == '<') {
        s.push(['m ',(r*(scr[1]-size)),', ',(r*(scr[2])),' l ',
        (r*(scr[1]+size/2)),', ',Math.round(r*(scr[2]-size*Math.sqrt(3)/2)),' l ',
        (r*(scr[1]+size/2)),', ',Math.round(r*(scr[2]+size*Math.sqrt(3)/2)),' x e '
        ].join(''));    
    }    
    return s;
}

JXG.VMLRenderer.prototype.updatePolygonePrimitive = function(node,el) {
    var minX = el.vertices[0].coords.scrCoords[1],
        maxX = el.vertices[0].coords.scrCoords[1],
        minY = el.vertices[0].coords.scrCoords[2],
        maxY = el.vertices[0].coords.scrCoords[2],
        i, 
        len = el.vertices.length,
        scr, x, y, 
        pStr = [];
        
    this.setAttr(node, 'stroked', 'false');
    for(i=1; i<len-1; i++) {
        scr = el.vertices[i].coords.scrCoords;
        if(scr[1] < minX) {
            minX = scr[1];
        }
        else if(scr[1] > maxX) {
            maxX = scr[1];
        }
        if(scr[2] < minY) {
            minY = scr[2];
        }
        else if(scr[2] > maxY) {
            maxY = scr[2];
        }
    }

    x = Math.round(maxX-minX); // Breite des umgebenden Rechtecks?
    y = Math.round(maxY-minY); // Hoehe des umgebenden Rechtecks?

    if (!isNaN(x) && !isNaN(y)) {
        node.style.width = x;
        node.style.height = y;
        this.setAttr(node, 'coordsize', x+','+y);
    }
     
    scr = el.vertices[0].coords.scrCoords;
    pStr.push(["m ",scr[1],",",scr[2]," l "].join(''));
    
    for(i=1; i<len-1; i++) {
        scr = el.vertices[i].coords.scrCoords;
        pStr.push(scr[1] + "," + scr[2]);
        if(i<len-2) {
            pStr.push(", ");
        }
    }
    pStr.push(" x e");

    this.setAttr(node, 'path',pStr.join(""));
};

JXG.VMLRenderer.prototype.appendChildPrimitive = function(node,level) {
    if (typeof level=='undefined') level = 0;   // For trace nodes    
    node.style.zIndex = level;
    /*
    switch (level) {
        case 'images': node.style.zIndex = "1"; break;
        case 'grid': node.style.zIndex = "1"; break;
        case 'angles': node.style.zIndex = "2"; break;
        case 'sectors': node.style.zIndex = "2"; break;
        case 'polygone': node.style.zIndex = "2"; break;
        case 'curves': node.style.zIndex = "4"; break; //2
        case 'circles': node.style.zIndex = "4"; break; //3
        case 'lines': node.style.zIndex = "4"; break;
        case 'arcs': node.style.zIndex = "4"; break;
        case 'points': node.style.zIndex = "5"; break;
    }
    */
    this.container.appendChild(node);
};

JXG.VMLRenderer.prototype.setPropertyPrimitive = function(node,key,val) {
    var keyVml = '', 
        node2, v;
        
    switch (key) {
        case 'stroke': keyVml = 'strokecolor'; break;
        case 'stroke-width': keyVml = 'strokeweight'; break;
        case 'stroke-dasharray': keyVml = 'dashstyle'; break;
    }
    if (keyVml!='') {
        v = this.eval(val);
        this.setAttr(node, keyVml, v);
    }
};

JXG.VMLRenderer.prototype.drawVerticalGrid = function(topLeft, bottomRight, gx, board) {
    var node = this.createPrimitive('path', 'gridx'),
        gridArr = [];
        
    while(topLeft.scrCoords[1] < bottomRight.scrCoords[1] + gx - 1) { 
        gridArr.push(' m ' + (this.resolution*topLeft.scrCoords[1]) + 
                     ', ' + 0 + 
                     ' l ' + (this.resolution*topLeft.scrCoords[1]) + 
                     ', ' + (this.resolution*board.canvasHeight)+' ');
        topLeft.setCoordinates(JXG.COORDS_BY_SCREEN, [topLeft.scrCoords[1] + gx, topLeft.scrCoords[2]]);   
    }
    this.updatePathPrimitive(node, gridArr, board);
    return node;
};

JXG.VMLRenderer.prototype.drawHorizontalGrid = function(topLeft, bottomRight, gy, board) {
    var node = this.createPrimitive('path', 'gridy'),
        gridArr = [];
    while(topLeft.scrCoords[2] <= bottomRight.scrCoords[2] + gy - 1) {
        gridArr.push(' m ' + 0 + 
                     ', ' + (this.resolution*topLeft.scrCoords[2]) + 
                     ' l ' + (this.resolution*board.canvasWidth) + 
                     ', ' + (this.resolution*topLeft.scrCoords[2])+' ');
        topLeft.setCoordinates(JXG.COORDS_BY_SCREEN, [topLeft.scrCoords[1], topLeft.scrCoords[2] + gy]);
    }
    this.updatePathPrimitive(node, gridArr, board);
    return node;
};

