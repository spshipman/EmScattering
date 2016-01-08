// Global namespace
var emScattering = emScattering || {};

// Use to exponentiate a diagonal matrix. Must pass in numeric.T type matrix.
emScattering.diagonalExp = function(matrix, power) {
    var diag = matrix.getDiag();
    var diagReal = diag.x, diagImag = diag.y, real, imag;
    real = numeric.mul(diagReal, power);
    imag = numeric.mul(diagImag, power);
    var zPrime = new numeric.T(real, imag);
    return numeric.T.diag(zPrime.exp());
};


/*----- PhotonicStructure1D -----*/
emScattering.PhotonicStructure1D = function(epsilon, mu, length) {
    this.epsilon = epsilon;
    this.mu = mu;
    this.length = length;
    this.N = length.length;
    this.overallTransferMatrix = [];
    this.transferMatrices = [];
    this.scatteringMatrix = [];
};

emScattering.PhotonicStructure1D.prototype.calcTransferMatrices = function(k0, kx, ky) {
    this.overallTransferMatrix = numeric.T.identity(4);
    this.transferMatrices = [];

    if (this.N == 1) {
        this.transferMatrices.push(this.overallTransferMatrix);
    }
    else {
        for (var i = 0; i < this.N-1; i++) {
            var eigenvectorsNextMatrix = this.eigenvectorsForLayer(i+1, k0, kx, ky);
            var eigenvectorsCurrentMatrix = this.eigenvectorsForLayer(i, k0, kx, ky);
            var eigenvalueCurrentMatrix = this.eigenvaluesForLayer(i, k0, kx, ky);

            var currentLayerLength;
            if (i == 0) {
                // The left ambient doesn't have a length as it isn't a true layer
                currentLayerLength = 0;
            }
            else {
                currentLayerLength = this.length[i];
            }

            var expEigenvalueMatrix = emScattering.diagonalExp(eigenvalueCurrentMatrix, currentLayerLength*k0);
            var currentLayerTransferMatrix = eigenvectorsNextMatrix.inv()
                .dot(eigenvectorsCurrentMatrix.dot(expEigenvalueMatrix));

            if (i == 0) {
                emScattering.print(expEigenvalueMatrix, "exp ");
                emScattering.print(eigenvectorsNextMatrix, "next ");
                emScattering.print(eigenvectorsCurrentMatrix, "current ");
            }

            this.transferMatrices.push(currentLayerTransferMatrix);
            this.overallTransferMatrix = currentLayerTransferMatrix.dot(this.overallTransferMatrix);
        }
    }
};

emScattering.PhotonicStructure1D.prototype.eigenvectorsForLayer = function(i, k0, kx, ky) {
    var eigenvectors = new numeric.T.rep([4,4], 0);
    var eigenvalueMatrix = this.eigenvaluesForLayer(i, k0, kx, ky);
    var forwardEigenvalue = eigenvalueMatrix.get([0,0]);

    var epsilon = this.epsilon[i];
    var mu = this.mu[i];

    var kx_norm = kx / k0;
    var ky_norm = ky / k0;

    // Convenience values for constructing the matrix
    var element1, element2, element3;
    element1 = new numeric.T(kx_norm*ky_norm/mu, 0);
    element1 = element1.div(forwardEigenvalue);

    element2 = new numeric.T((mu*epsilon - Math.pow(kx_norm, 2))/mu,0);
    element2 = element2.div(forwardEigenvalue);

    element3 = new numeric.T((mu*epsilon - Math.pow(ky_norm, 2))/mu,0);
    element3 = element3.div(forwardEigenvalue);

    // Eigenvector 1 - forward prop
    eigenvectors.set([0, 0], new numeric.T(1, 0));
    eigenvectors.set([2, 0], element1);
    eigenvectors.set([3, 0], element3.dot(-1));

    // Eigenvector 2 - forward prop
    eigenvectors.set([1, 1], new numeric.T(1, 0));
    eigenvectors.set([2, 1], element2);
    eigenvectors.set([3, 1], element1.dot(-1));

    // Eigenvector 3 - backward prop
    eigenvectors.set([0, 2], new numeric.T(1, 0));
    eigenvectors.set([2, 2], element1.dot(-1));
    eigenvectors.set([3, 2], element3);

    // Eigenvector 4 - backward prop
    eigenvectors.set([1, 3], new numeric.T(1, 0));
    eigenvectors.set([2, 3], element2.dot(-1));
    eigenvectors.set([3, 3], element1);

    return eigenvectors;
};

emScattering.PhotonicStructure1D.prototype.eigenvaluesForLayer = function(i, k0, kx, ky) {
    var epsilon = this.epsilon[i];
    var mu = this.mu[i];
    var kxNormSquared = Math.pow(kx / k0, 2);
    var kyNormSquared = Math.pow(ky / k0, 2);

    var n = mu * epsilon;
    var eigenvalueSquared = (kxNormSquared * kyNormSquared - (n - kxNormSquared) * (n - kyNormSquared)) / n;

    var forwardEigenvalue, backwardEigenvalue;
    if (eigenvalueSquared >= 0) {
        forwardEigenvalue = new numeric.T(Math.sqrt(eigenvalueSquared), 0);
    }
    else {
        forwardEigenvalue = new numeric.T(0, Math.sqrt(-1*eigenvalueSquared));
    }
    backwardEigenvalue = forwardEigenvalue.dot(-1);

    var eigenvalues = numeric.T.rep([4,4], 0);
    eigenvalues.set([0,0], forwardEigenvalue);
    eigenvalues.set([1,1], forwardEigenvalue);
    eigenvalues.set([2,2], backwardEigenvalue);
    eigenvalues.set([3,3], backwardEigenvalue);
    return eigenvalues;
};

emScattering.PhotonicStructure1D.prototype.calcScatteringMatrix = function(k0, kx, ky) {
    // Initialize global scattering matrix blocks -- perfect transmission and no reflection
    var globalScatteringMatrix = new emScattering.ScatteringMatrix(new numeric.T.rep([2, 2], 0), new numeric.T.identity(2),
        new numeric.T.identity(2), new numeric.T.rep([2, 2], 0));

    // Layer scattering matrices
    for (var i = 1; i < this.N-1; i++) {
        var currentLayerScatteringMatrix = emScattering.ScatteringMatrix.constructScatteringMatrix(this, i, k0, kx, ky);
        globalScatteringMatrix = globalScatteringMatrix.combine(currentLayerScatteringMatrix);
    }

    // Connect to reflection region - left side
    var reflectionSideScatteringMatrix = emScattering.ScatteringMatrix.constructScatteringMatrix(this, 0, k0, kx, ky);
    globalScatteringMatrix = reflectionSideScatteringMatrix.combine(globalScatteringMatrix);

    // Connect to transmission region - right side
    var transmissionSideScatteringMatrix = emScattering.ScatteringMatrix.constructScatteringMatrix(this, this.N-1, k0, kx, ky);
    globalScatteringMatrix = globalScatteringMatrix.combine(transmissionSideScatteringMatrix);

    this.scatteringMatrix = globalScatteringMatrix;
    return this.scatteringMatrix;
};

emScattering.PhotonicStructure1D.prototype.getTransmittance = function(k0 ,kx, ky) {
    var kxNorm = kx / k0;
    var kyNorm = ky / k0;

    var kzTransmitted = this.getKz(this.N-1, k0, kx, ky);
    var kzReflected = this.getKz(0, k0, kx, ky);

    this.calcScatteringMatrix(k0, kx, ky);

    var incoming = new numeric.T.rep([2, 1], 0);
    var n = new numeric.T(1/Math.sqrt(2), 0);
    incoming.set([0, 0], n);
    incoming.set([1, 0], n);
    var incomingEz = (n.mul(kxNorm).add(n.mul(kyNorm))).mul(-1).div(kzReflected);
    var incomingEnergySquared = 2*Math.pow(incoming.norm2(), 2) + Math.pow(incomingEz.norm2(), 2);

    var s21 = this.scatteringMatrix.s21;

    var transmitted = s21.dot(incoming);
    emScattering.print(transmitted, "trn");
    var transmittedEx = transmitted.get([0, 0]);
    var transmittedEy = transmitted.get([1, 0]);
    var transmittedEz = transmittedEx.mul(kxNorm).add(transmittedEy.mul(kyNorm)).mul(-1).div(kzTransmitted);
    var transmittedEnergySquared =  Math.pow(transmittedEx.norm2(), 2) +  Math.pow(transmittedEy.norm2(), 2) +  Math.pow(transmittedEz.norm2(), 2);

    var tMult = kzTransmitted.div(kzReflected).mul(this.mu[0]/this.mu[this.N-1]);
    return transmittedEnergySquared / incomingEnergySquared * tMult.x;
};

emScattering.PhotonicStructure1D.prototype.getReflectance = function(k0, kx, ky) {
    var kxNorm = kx / k0;
    var kyNorm = ky / k0;

    var kzReflected = this.getKz(0, k0, kx, ky);

    this.calcScatteringMatrix(k0, kx, ky);

    var incoming = new numeric.T.rep([2, 1], 0);
    var n = new numeric.T(1/Math.sqrt(2), 0);
    incoming.set([0, 0], n);
    incoming.set([1, 0], n);
    var incomingEz = (n.mul(kxNorm).add(n.mul(kyNorm))).mul(-1).div(kzReflected);
    var incomingEnergySquared = Math.pow(incoming.norm2(), 2) + Math.pow(incomingEz.norm2(), 2);

    var s11 = this.scatteringMatrix.s11;

    var reflected = s11.dot(incoming);
    emScattering.print(reflected);
    var reflectedEx = reflected.get([0, 0]);
    var reflectedEy = reflected.get([1, 0]);
    var reflectedEz = reflectedEx.mul(kxNorm).add(reflectedEy.mul(kyNorm)).mul(-1).div(kzReflected);
    var reflectedEnergySquared = Math.pow(reflectedEx.norm2(), 2) + Math.pow(reflectedEy.norm2(), 2) + Math.pow(reflectedEz.norm2(), 2);

    return reflectedEnergySquared / incomingEnergySquared;
};

emScattering.PhotonicStructure1D.prototype.polarizationVector = function(k0, kx, ky) {
    var kIncoming = numeric.T.rep([3, 1], 0);
    kIncoming.set([0, 0], kx);
    kIncoming.set([0, 1], ky);
    var kz = this.getKz(0, k0, kx, ky).mul(k0);
    kIncoming.set([0, 2], kz);

    var surfaceNormal = numeric.T.rep([3, 1], 0);
    surfaceNormal.set([2, 0], 1);

    var tePolarizationVector = new numeric.T.rep([3, 1], 0), tmPolarizationVector = new numeric.T.rep([3, 1], 0);
    if (kx == 0 && ky == 0) {
        tePolarizationVector.set([1, 0], 1);
        tmPolarizationVector.set([0, 0], kz);
        tmPolarizationVector.set([1, 0], 0);
        tmPolarizationVector.set([2, 0], -1*kx);
        emScattering.print(tmPolarizationVector, "tm vec");
        tmPolarizationVector = tmPolarizationVector.mul(1/tmPolarizationVector.norm2());
    }
    else {
        tePolarizationVector.set([0, 0], kz);
        tePolarizationVector.set([2, 0], -1*kx);
        tePolarizationVector = tePolarizationVector.div(tePolarizationVector.norm2());

        var ax = tePolarizationVector.get([0, 0]);
        var az = tePolarizationVector.get([2, 0]);
        tmPolarizationVector.set([0, 0], az.mul(-1*ky));
        tmPolarizationVector.set([1, 0], az.mul(kx).sub(ax.mul(kz)));
        tmPolarizationVector.set([2, 0], ax.mul(ky));
        tmPolarizationVector = tmPolarizationVector.div(tmPolarizationVector.norm2());
    }
    emScattering.print(tePolarizationVector, "te ");
    emScattering.print(tmPolarizationVector.norm2(), "tm ");

    var polarization = tePolarizationVector.add(tmPolarizationVector);
    polarization = polarization.div(polarization.norm2());
    return polarization;
}

// Normalized kz
emScattering.PhotonicStructure1D.prototype.getKz = function(i, k0, kx, ky) {
    var kxNorm = kx / k0;
    var kyNorm = ky / k0;
    var kzSquared = this.mu[i] * this.epsilon[i] - Math.pow(kxNorm, 2) - Math.pow(kyNorm, 2);

    var kz;
    if (kzSquared >= 0) {
        kz = new numeric.T(Math.sqrt(kzSquared), 0);
    }
    else {
        kz = new numeric.T(0, Math.sqrt(-1*kzSquared));
    }
    return kz;
};

// TODO: Bug in transfer matrix calculation (or in scattering matrix) as the two don't predict the same U value
emScattering.PhotonicStructure1D.prototype.determineField = function(k0, kx, ky, J) {
    var numPoints = 100;    // Number of points calculated per unit length
    var _Ex = [], _Ey = [], _Hx = [], _Hy = [], _z = [];

    this.calcTransferMatrices(k0, kx, ky);

    var U = this.calcScatteringMatrix(k0, kx, ky).dot(J);

    var cx = [J[0], J[1], U.x[0], U.x[1]];
    var cy = [0, 0, U.y[0], U.y[1]];
    var c = new numeric.T(cx, cy);   // Coefficients in the current layer

    // Determine the fields in each layer
    var eigenvectors, eigenvalues, layerNormZ, result, i, j, expDiag;
    var currentLeftZ = 0, currentRightZ = this.length[0];
    for (i = 0; i < this.N; i++) {
        if (i == 0) {
            layerNormZ= numeric.linspace(-1*k0*this.length[i], 0,
                this.length[i]*numPoints);
        }
        else {
            layerNormZ = numeric.linspace(0, k0*this.length[i],
                this.length[i]*numPoints);
        }
        eigenvectors = this.eigenvectorsForLayer(i, k0, kx, ky);
        eigenvalues = this.eigenvaluesForLayer(i, k0, kx, ky);

        for (j = 0, M = layerNormZ.length; j < M; j++) {
            expDiag = emScattering.diagonalExp(eigenvalues, layerNormZ[j]);
            result = eigenvectors.dot(expDiag.dot(c));
            _Ex.push(result.x[0]);
            _Ey.push(result.x[1]);
            _Hx.push(result.x[2]);
            _Hy.push(result.x[3]);
        }

        _z = _z.concat(numeric.linspace(currentLeftZ, currentRightZ,
            numPoints*this.length[i]));
        if (i+1 < this.N) {
            currentLeftZ += this.length[i];
            currentRightZ += this.length[i+1];
            c = this.transferMatrices[i].dot(c);
        }
    }
    //emScattering.print(U, "U from S: ");
    //emScattering.print(c, "Outgoing on left (last 2 elements): ");
    return {z: _z, Ex: _Ex, Ey: _Ey, Hx: _Hx, Hy: _Hy};
};


emScattering.ScatteringMatrix = function(s11, s12, s21, s22) {
    this.s11 = s11;
    this.s12 = s12;
    this.s21 = s21;
    this.s22 = s22;
};

// Construct the scattering matrix for the given layer i surrounded by a gap layer
emScattering.ScatteringMatrix.constructScatteringMatrix = function(structure, i, k0, kx, ky) {
    var kxNorm = kx / k0;
    var kyNorm = ky / k0;
    var kz = structure.getKz(i, k0, kx, ky);

    // Gap layer parameters that is chosen so that any wave propagates within it. It's vacuum.
    var gapLayerQMatrix = numeric.T.rep([2, 2], 0);
    gapLayerQMatrix.set([0, 0], new numeric.T(kxNorm*kyNorm, 0));
    gapLayerQMatrix.set([0, 1], new numeric.T(1+Math.pow(kyNorm,2)));
    gapLayerQMatrix.set([1, 0], new numeric.T(-1*Math.pow(kxNorm,2)-1));
    gapLayerQMatrix.set([1, 1], new numeric.T(-1*kxNorm*kyNorm));

    var gapLayerVMatrix = gapLayerQMatrix.mul(new numeric.T(0, -1));

    // Layer parameters
    var mu = structure.mu[i];
    var epsilon = structure.epsilon[i];

    var qMatrix = numeric.T.rep([2, 2], 0);
    qMatrix.set([0, 0], new numeric.T(kxNorm*kyNorm/mu, 0));
    qMatrix.set([0, 1], new numeric.T((mu*epsilon-Math.pow(kxNorm,2))/mu,0));
    qMatrix.set([1, 0], new numeric.T((Math.pow(kyNorm,2)-mu*epsilon)/mu,0));
    qMatrix.set([1, 1], new numeric.T((-1*kxNorm*kyNorm)/mu,0));

    var omegaMatrix = numeric.T.identity(2).mul(new numeric.T(0, 1));
    omegaMatrix.mul(kz);

    var vMatrix = qMatrix.dot(omegaMatrix.inv());

    // Scattering matrix for layer i
    var aMatrix = numeric.T.identity(2).add(vMatrix.inv().dot(gapLayerVMatrix));
    var bMatrix = numeric.T.identity(2).sub(vMatrix.inv().dot(gapLayerVMatrix));

    var eigenvaluesMatrix = structure.eigenvaluesForLayer(i, k0, kx, ky);
    var xMatrix = numeric.T.rep([2, 2], 0);
    xMatrix.set([0, 0], eigenvaluesMatrix.get([0, 0]));
    xMatrix.set([1, 1], eigenvaluesMatrix.get([0, 0]));
    xMatrix = emScattering.diagonalExp(xMatrix, k0*structure.length[i]);

    var dMatrix = aMatrix.sub(xMatrix.dot(bMatrix).dot(aMatrix.inv()).dot(xMatrix).dot(bMatrix));

    var s11 = dMatrix.inv().dot(xMatrix.dot(bMatrix).dot(aMatrix.inv()).dot(xMatrix).dot(aMatrix).sub(bMatrix));
    var s22 = s11;
    var s12 = dMatrix.inv().dot(xMatrix).dot(aMatrix.sub(bMatrix.dot(aMatrix.inv()).dot(bMatrix)));
    var s21 = s12;
    return new emScattering.ScatteringMatrix(s11, s12, s21, s22);
};

// Scattering matrix for AB (star product of A * B) where A is 'this'
emScattering.ScatteringMatrix.prototype.combine = function(B) {
    var D = this.s12.dot((numeric.T.identity(2).sub(B.s11.dot(this.s22))).inv());
    var F = B.s21.dot((numeric.T.identity(2).sub(this.s22.dot(B.s11))).inv());

    var s11 = this.s11.add(D.dot(B.s11).dot(this.s21));
    var s12 = D.dot(B.s12);
    var s21 = F.dot(this.s21);
    var s22 = B.s22.add(F.dot(this.s22).dot(B.s12));
    return new emScattering.ScatteringMatrix(s11, s12, s21, s22);
};

emScattering.ScatteringMatrix.prototype.dot = function(B) {
    var S = numeric.T.rep([4,4], 0);
    S.setBlock([2, 0], [3, 1], this.s21);
    S.setBlock([2, 2], [3, 3], this.s22);
    S.setBlock([0, 0], [1, 1], this.s11);
    S.setBlock([0, 2], [1, 3], this.s12);
    return S.dot(B);
}


emScattering.print = function(mat, pre) {
    if (pre == null)
        pre = "";
    //document.write(pre + numeric.prettyPrint(mat) + "<br>");
};

emScattering.printFields = function(fields) {
    document.write("<table>");
    document.write("<tr><td>z</td> <td>Ex</td> <td>Ey</td> <td>Hx</td> <td>Hy</td></tr>");
    for (var i = 0, N = fields.z.length; i < N; i++) {
        document.write("<tr>");
        document.write("<td>" + fields.z[i] + "</td>");
        document.write("<td>" + fields.Ex[i] + "</td>");
        document.write("<td>" + fields.Ey[i] + "</td>");
        document.write("<td>" + fields.Hx[i] + "</td>");
        document.write("<td>" + fields.Hy[i] + "</td>");
        document.write("</tr>");
    }
    document.write("</table>");
};

var k0 = 3;
var kx = 0;
var ky = 0;
var epsilon = [1, 5, 1];
var mu = [1, 5, 1];
var length = [0, 5, 1];
var structure = new emScattering.PhotonicStructure1D(epsilon, mu, length);

//emScattering.print(structure.eigenvaluesForLayer(1,k0,kx,ky));
//emScattering.print(structure.eigenvectorsForLayer(1,k0,kx,ky));
//
//structure.calcTransferMatrices(k0, kx, ky);
//emScattering.print(structure.transferMatrices[1]);
//
//var s = structure.calcScatteringMatrix(k0, kx, ky);
//emScattering.print(s, "s ");
//emScattering.print(s.dot(s.transjugate()), "eye(4) ");
//
emScattering.print(structure.getTransmittance(k0,kx,ky)+structure.getReflectance(k0,kx,ky),"T ");
emScattering.print(structure.polarizationVector(k0, kx, ky));
emScattering.printFields(structure.determineField(k0, kx, ky, [1,0,0,0]));
