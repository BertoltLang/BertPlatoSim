
#include "Camera.h"



pair<double, double> Camera::newSkyToFocalPlaneCoordinates(double raStar, double decStar, double startTime, double time, bool useJitter, bool useDrift)
{
    // create the initial spacecraft coordinate system conversion matrix according 
    // to PLATO-DLR-PL-TN-016_i1_2draft by Denis Griessbach

    // the matrix transformations consist of rotations as well as translations, 
    // which is why 4x4 matrices are used

    // S/C reference frame

    // the z-axis is the pointing coordinate of the spacecraft, the y-axis is lying in the equatorial plane 
    // and the x-axis is calculated
    arma::vec zSC = {cos(decPlatform) * cos(alphaPlatform), cos(decPlatform) * sin(alphaPlatform), sin(decPlatform)};
    arma::vec ySC = {cos(alphaPlatform + PI/2), sin(alphaPlatform + PI/2), 0.0};
    arma::vec xSC = cross(ySC, zSC);

    arma::mat rotSCtoEQ;
    rotSCtoEQ       <<  xSC[0] << ySC[0] << zSC[0] << 0.0 << arma::endr
                    <<  xSC[1] << ySC[1] << zSC[1] << 0.0 << arma::endr
                    <<  xSC[2] << ySC[2] << zSC[2] << 0.0 << arma::endr
                    <<    0.0  <<   0.0  <<   0.0  << 1.0 << arma::endr;

    // include a rotation around the z - axis to align the x - axis to the average direction of the sun

    const double sunDirectionAngle = deg2rad(-9.047422535);

    //create the rotation matrix
    double sinSun = sin(sunDirectionAngle);
    double cosSun = cos(sunDirectionAngle);

    arma::mat rotSun;
    rotSun          <<   cosSun << -sinSun << 0.0 << 0.0 << arma::endr
                    <<   sinSun <<  cosSun << 0.0 << 0.0 << arma::endr
                    <<    0.0   <<   0.0   << 1.0 << 0.0 << arma::endr
                    <<    0.0   <<   0.0   << 0.0 << 1.0 << arma::endr;

    arma::mat rotEQtoSC;

    // include jitter, if it's considered

    if (useJitter)
    {

        //include jitter
        double yaw=0.0, pitch=0.0, roll=0.0;
        double timeInterval = time - startTime;
        tie(yaw, pitch, roll) = telescope.getNextYawPitchRoll(timeInterval);

        // Some handy abbreviations

        const double sinYaw   = sin(yaw);
        const double cosYaw   = cos(yaw);
        const double sinPitch = sin(pitch);
        const double cosPitch = cos(pitch);
        const double sinRoll  = sin(roll);
        const double cosRoll  = cos(roll);

        // The rotation matrices

        arma::mat rotYaw;
        rotYaw      << 1.0  <<    0.0   <<     0.0   << 0.0 << arma::endr
                    << 0.0  <<  cosYaw  <<  -sinYaw  << 0.0 << arma::endr
                    << 0.0  <<  sinYaw  <<   cosYaw  << 0.0 << arma::endr
                    << 0.0  <<    0.0   <<     0.0   << 1.0 << arma::endr;

        arma::mat rotPitch;
        rotPitch    <<  cosPitch  <<  0.0  <<  sinPitch  << 0.0 << arma::endr
                    <<   0.0      <<  1.0  <<    0.0     << 0.0 << arma::endr
                    << -sinPitch  <<  0.0  <<  cosPitch  << 0.0 << arma::endr
                    <<   0.0      <<  0.0  <<    0.0     << 1.0 << arma::endr;

        arma::mat rotRoll;
        rotRoll     << cosRoll  <<  -sinRoll  <<  0.0  << 0.0 << arma::endr
                    << sinRoll  <<   cosRoll  <<  0.0  << 0.0 << arma::endr
                    <<  0.0     <<     0.0    <<  1.0  << 0.0 << arma::endr
                    <<  0.0     <<     0.0    <<  0.0  << 1.0 << arma::endr;

        //rotEQtoSC = rotRoll.t() * rotPitch.t() * rotYaw.t() * rotSCtoEQ.t();

        rotEQtoSC = rotSCtoEQ.t(); 

        rotEQtoSC = rotRoll * rotPitch * rotYaw * rotSun * rotEQtoSC;

    }
    else
    {
        rotEQtoSC = rotSun * rotSCtoEQ.t();
    }

    // the payload module reference frame

    // the PLM frame coincides with the SC frame, except for a small mounting  offset, 
    // which is a translation on the z-axis in positive direction. This is just for completeness, since
    // an offset in z - direction doesn't influence the picture generation

    arma::vec distSCtoPLM = {0.0, 0.0, 0.0};

    arma::mat rotSCtoPLM;
    rotSCtoPLM      << 1.0 << 0.0 << 0.0 << distSCtoPLM[0] << arma::endr
                    << 0.0 << 1.0 << 0.0 << distSCtoPLM[1] << arma::endr
                    << 0.0 << 0.0 << 1.0 << distSCtoPLM[2] << arma::endr
                    << 0.0 << 0.0 << 0.0 <<      1.0       << arma::endr;

    // include the thermoelastic drift TED, if it's considered

    if (useDrift)
    {
        //include drift
        double yaw=0.0, pitch=0.0, roll=0.0;
        double timeInterval = time - internalTime;
        tie(yaw, pitch, roll) = telescope.getThermalDrift(timeInterval);

        // Some handy abbreviations

        const double sinYaw   = sin(yaw);
        const double cosYaw   = cos(yaw);
        const double sinPitch = sin(pitch);
        const double cosPitch = cos(pitch);
        const double sinRoll  = sin(roll);
        const double cosRoll  = cos(roll);

        // The rotation matrices

        arma::mat rotYaw;
        rotYaw      << 1.0  <<    0.0   <<     0.0   << 0.0 << arma::endr
                    << 0.0  <<  cosYaw  <<  -sinYaw  << 0.0 << arma::endr
                    << 0.0  <<  sinYaw  <<   cosYaw  << 0.0 << arma::endr
                    << 0.0  <<    0.0   <<     0.0   << 1.0 << arma::endr;

        arma::mat rotPitch;
        rotPitch    <<  cosPitch  <<  0.0  <<  sinPitch  << 0.0 << arma::endr
                    <<   0.0      <<  1.0  <<    0.0     << 0.0 << arma::endr
                    << -sinPitch  <<  0.0  <<  cosPitch  << 0.0 << arma::endr
                    <<   0.0      <<  0.0  <<    0.0     << 1.0 << arma::endr;

        arma::mat rotRoll;
        rotRoll     << cosRoll  <<  -sinRoll  <<  0.0  << 0.0 << arma::endr
                    << sinRoll  <<   cosRoll  <<  0.0  << 0.0 << arma::endr
                    <<  0.0     <<     0.0    <<  1.0  << 0.0 << arma::endr
                    <<  0.0     <<     0.0    <<  0.0  << 1.0 << arma::endr;


        rotSCtoPLM = rotRoll * rotPitch * rotYaw * rotSCtoPLM;
    }


    // the camera/FPA reference frame

    // the x - axis should point towards the orientation screw, which can be achieved 
    // by a rotation around the z  -axis

    const double screwAngle = 0.0;

    double sinScrew = sin(screwAngle);
    double cosScrew = cos(screwAngle);

    arma::mat rotScrew;
    rotScrew        <<   cosScrew << -sinScrew << 0.0 << 0.0 << arma::endr
                    <<   sinScrew <<  cosScrew << 0.0 << 0.0 << arma::endr
                    <<      0.0   <<     0.0   << 1.0 << 0.0 << arma::endr
                    <<      0.0   <<     0.0   << 0.0 << 1.0 << arma::endr;     

    // the z - axis points towards the line of sight which is achieved by a 
    // rotation of the z axis around the x - axis over the tilt - angle 

    double sinTilt = sin(tiltAngle);
    double cosTilt = cos(tiltAngle);

    arma::mat rotTilt;
    rotTilt         << 1.0 << 0.0     <<    0.0   << 0.0 << arma::endr
                    << 0.0 << cosTilt << -sinTilt << 0.0  << arma::endr
                    << 0.0 << sinTilt <<  cosTilt << 0.0 << arma::endr
                    << 0.0 <<   0.0   <<    0.0   << 1.0 << arma::endr;  

    // the next step is the rotation over the focal plane angle

    double sinOrient = sin(focalPlaneAngle);
    double cosOrient = cos(focalPlaneAngle);

    arma::mat rotOrient;
    rotOrient       <<   cosOrient <<  -sinOrient << 0.0 << 0.0 << arma::endr
                    <<   sinOrient <<   cosOrient << 0.0 << 0.0 << arma::endr
                    <<     0.0     <<     0.0     << 1.0 << 0.0 << arma::endr
                    <<     0.0     <<     0.0     << 0.0 << 1.0 << arma::endr; 

    //rotEQtoSC = rotRoll * rotPitch * rotYaw * rotOrient * rotEQtoSC;                


    arma::mat rotPLMtoCAM = rotOrient * rotTilt * rotScrew;


    // combine all transformation matrices

    arma::mat rotEQtoCAM = rotPLMtoCAM * rotSCtoPLM * rotEQtoSC;

    // create the cartesian star coordinates
    arma::vec starEQ = {cos(decStar) * cos(raStar), cos(decStar) * sin(raStar), sin(decStar), 1};

    // calculation of the star coordinates in the CAM reference frame   




    //printf("mit Jitter: \n");
    //cout << rotEQtoCAM;
    //rotEQtoCAM.print();


    arma::vec starCAM = rotEQtoCAM * starEQ;

    //calculate the mm values of the coordinates and return them

    const double xFPmm = - focalLength * starCAM(0)/starCAM(2);
    const double yFPmm = - focalLength * starCAM(1)/starCAM(2);

    return make_pair(xFPmm, yFPmm);

}


//newFocalPlaneToSkyCoordinates replaces focalPlaneToSkyCoordinates

pair<double, double> Camera::newFocalPlaneToSkyCoordinates(double xFPprime, double yFPprime)
{    
    arma::vec zSC = {cos(decPlatform) * cos(alphaPlatform), cos(decPlatform) * sin(alphaPlatform), sin(decPlatform)};
    arma::vec ySC = {cos(alphaPlatform + PI/2), sin(alphaPlatform + PI/2), 0.0};
    arma::vec xSC = cross(ySC, zSC);

    arma::mat rotSCtoEQ;
    rotSCtoEQ       <<  xSC[0] << ySC[0] << zSC[0] << 0.0 << arma::endr
                    <<  xSC[1] << ySC[1] << zSC[1] << 0.0 << arma::endr
                    <<  xSC[2] << ySC[2] << zSC[2] << 0.0 << arma::endr
                    <<    0.0  <<   0.0  <<   0.0  << 1.0 << arma::endr;

    // include a rotation around the z - axis to align the x - axis to the average direction of the sun

    const double sunDirectionAngle = deg2rad(-9.047422535);

    //create the rotation matrix
    double sinSun = sin(sunDirectionAngle);
    double cosSun = cos(sunDirectionAngle);

    arma::mat rotSun;
    rotSun          <<   cosSun << -sinSun << 0.0 << 0.0 << arma::endr
                    <<   sinSun <<  cosSun << 0.0 << 0.0 << arma::endr
                    <<    0.0   <<   0.0   << 1.0 << 0.0 << arma::endr
                    <<    0.0   <<   0.0   << 0.0 << 1.0 << arma::endr;

    rotSCtoEQ = rotSCtoEQ * rotSun.t();


    // the payload module reference frame

    // the PLM frame coincides with the SC frame, except for a small mounting  offset, 
    // which is a translation on the z-axis in positive direction. This is just for completeness, since
    // an offset in z - direction doesn't influence the picture generation

    arma::vec distSCtoPLM = {0.0, 0.0, 0.0};

    arma::mat rotPLMtoSC;
    rotPLMtoSC      << 1.0 << 0.0 << 0.0 << -distSCtoPLM[0] << arma::endr
                    << 0.0 << 1.0 << 0.0 << -distSCtoPLM[1] << arma::endr
                    << 0.0 << 0.0 << 1.0 << -distSCtoPLM[2] << arma::endr
                    << 0.0 << 0.0 << 0.0 <<      1.0       << arma::endr;


    // the camera/FPA reference frame

    // the x - axis should point towards the orientation screw, which can be achieved 
    // by a rotation around the z  -axis

    const double screwAngle = 0.0;

    double sinScrew = sin(screwAngle);
    double cosScrew = cos(screwAngle);

    arma::mat rotScrew;
    rotScrew        <<   cosScrew << -sinScrew << 0.0 << 0.0 << arma::endr
                    <<   sinScrew <<  cosScrew << 0.0 << 0.0 << arma::endr
                    <<      0.0   <<     0.0   << 1.0 << 0.0 << arma::endr
                    <<      0.0   <<     0.0   << 0.0 << 1.0 << arma::endr;     

    // the z - axis points towards the line of sight which is achieved by a 
    // rotation of the z axis around the x - axis over the tilt - angle 

    double sinTilt = sin(tiltAngle);
    double cosTilt = cos(tiltAngle);

    arma::mat rotTilt;
    rotTilt         << 1.0 << 0.0     <<    0.0   << 0.0 << arma::endr
                    << 0.0 << cosTilt << -sinTilt << 0.0  << arma::endr
                    << 0.0 << sinTilt <<  cosTilt << 0.0 << arma::endr
                    << 0.0 <<   0.0   <<    0.0   << 1.0 << arma::endr;  

    // the next step is the rotation over the focal plane angle

    double sinOrient = sin(focalPlaneAngle);
    double cosOrient = cos(focalPlaneAngle);

    arma::mat rotOrient;
    rotOrient       <<   cosOrient <<  -sinOrient << 0.0 << 0.0 << arma::endr
                    <<   sinOrient <<   cosOrient << 0.0 << 0.0 << arma::endr
                    <<     0.0     <<     0.0     << 1.0 << 0.0 << arma::endr
                    <<     0.0     <<     0.0     << 0.0 << 1.0 << arma::endr; 


    arma::mat rotCAMtoPLM = rotScrew.t() * rotTilt.t() * rotOrient.t();

    // combine all transformation matrices

    arma::mat rotCAMtoEQ = rotSCtoEQ * rotPLMtoSC * rotCAMtoPLM;

    //reverse the projection
    arma::vec starCAM = {xFPprime/-focalLength, yFPprime/-focalLength, 1.0, 1.0};

    // convert the CAM star coordinates to EQ star coordinates
    arma::vec starEQ = rotCAMtoEQ * starCAM;

    //convert the cartesian coordinates to equatorial coordinates
    double norm = sqrt(starEQ(0)*starEQ(0) + starEQ(1)*starEQ(1) + starEQ(2)*starEQ(2)); 
    double decStar = Constants::PI/2.0 - acos(starEQ(2)/norm);
    double raStar = atan2(starEQ(1), starEQ(0));

    if (raStar < 0.0)
    {
        raStar += 2.*Constants::PI;
    }

    // That's it!

    return make_pair(raStar, decStar);

}







/**
 * \brief      Constructor
 *
 * \param configParam   Configuration parameters for the Camera
 * \param hdf5file      HFD5 file to write the camera information to
 * \param telescope     Telescope on which the camera is mounted on
 * \param sky           Sky object to query which stars we are seeing
 */

Camera::Camera(ConfigurationParameters &configParam, HDF5File &hdf5file, Telescope &telescope, Sky &sky)
: HDF5Writer(hdf5file), telescope(telescope), sky(sky), internalTime(0.0), includeFieldDistortion(true), 
  userGivenSkyBackground(-1.0), fluxOfV0Star(0.0)
{
    // Create the groups in the HDF5 file where the different arrays generated by Camera will be saved.

    initHDF5Groups();

    // Parse the parameters from the configuration file.

    configure(configParam);

    // Initialize and load the PSF. This will open the PSF HDF5 file and perform some basic checking, 
    // but the psf itself will only be loaded by the psf.select(radius) method.

    psf = new PointSpreadFunction(configParam, hdf5file);

    // Initialize the polynomials that describes the field distortion of the camera.

    polynomial = Polynomial1D(polynomialDegree, polynomialCoefficients);
    inversePolynomial = Polynomial1D(polynomialDegree, inversePolynomialCoefficients);
}









/**
 * \brief  Destructor, free memory for psf
 */

Camera::~Camera()
{
    flushOutput();
    delete psf;
}













/**
 * \brief Creates the group(s) in the HDF5 file where the star information will be stored. 
 *        These group(s) have to be created once, at the very beginning.
 */

void Camera::initHDF5Groups()
{
    Log.debug("Camera: initialising HDF5 groups");

    hdf5File.createGroup("/StarPositions");
    hdf5File.createGroup("/Background");
}








/**
 * \brief      Collect and return the IDs of all stars that fall within the subField
 * 
 * \details    Note that this method pulls the detected stars from a map that is filled by 
 *             exposeDetector for each exposure. So, depending on when this method is called, 
 *             the returned set might be empty or incomplete.
 *             
 * \return     a set of unique star IDs that were detected in the subField
 */
set<unsigned int> Camera::getAllStarIDs()
{
    set<unsigned int> allStarIDs;            // A set<> stores only unique members

    for(auto timeMapPair: detectedStarInfo)
    {
        for (auto idArrayPair: timeMapPair.second)
        {
            allStarIDs.insert(idArrayPair.first);
        }
    }

    return allStarIDs;

}










/**
 * \brief Write all recorded information to the HDF5 output file
 * 
 */

void Camera::flushOutput()
{
    Log.info("Camera: Flushing output to HDf5 file.");

    // Extract and save the time points of all exposures 
    // Note: keyValuePair is (key, value) pair, where key is also a pair consisting of the startTime and StarID

    vector<double> time;
    for(auto keyValuePair: detectedStarInfo) time.push_back(keyValuePair.first);
    if (!time.empty())
    {
        hdf5File.writeArray("StarPositions/", "Time", time.data(), time.size());
    }
    else
    {
        Log.warning("Camera: No star positions to write to HDF5 file.");
    }


    // For each of the exposures, make a subgroup and write the position and flux of all detected stars.
    // Because some stars at the edge may jitter in and out of the subfield from one exposure to the other,
    // the written arrays may not be equally long for each exposure.


    for (int n = 0; n < time.size(); n++)
    {
        // Make the subgroup group

        stringstream myStream;
        myStream << "Exposure" << setfill('0') << setw(6) << n;
        const string exposureGroupName = "/StarPositions/" + myStream.str();
        hdf5File.createGroup(exposureGroupName);

        // Collect the different time series. For the positions, we only compute the sum, so we still need
        // to divide by N to compute the average, where N is the number of times the star was detected to be
        // in the subfield during an exposure.

        vector<unsigned int> starIDs;
        vector<double> xFPmm;
        vector<double> yFPmm;
        vector<double> rowPix;
        vector<double> colPix;
        vector<double> flux;

        for(auto keyValuePair: detectedStarInfo[time[n]])
        {
            const unsigned int starID = keyValuePair.first;
            starIDs.push_back(starID);                       // list of starIDs for this exposure only
            xFPmm.push_back(detectedStarInfo[time[n]][starID][0] / detectedStarInfo[time[n]][starID][5]);
            yFPmm.push_back(detectedStarInfo[time[n]][starID][1] / detectedStarInfo[time[n]][starID][5]);
            rowPix.push_back(detectedStarInfo[time[n]][starID][2] / detectedStarInfo[time[n]][starID][5]);
            colPix.push_back(detectedStarInfo[time[n]][starID][3] / detectedStarInfo[time[n]][starID][5]);
            flux.push_back(detectedStarInfo[time[n]][starID][4]);
        }

        // Write the time series to HDF5

        if(!starIDs.empty())
        {
            hdf5File.writeArray(exposureGroupName, "starID", starIDs.data(), starIDs.size());
            hdf5File.writeArray(exposureGroupName, "xFPmm",  xFPmm.data(),   xFPmm.size());
            hdf5File.writeArray(exposureGroupName, "yFPmm",  yFPmm.data(),   yFPmm.size());
            hdf5File.writeArray(exposureGroupName, "rowPix", rowPix.data(),  rowPix.size());
            hdf5File.writeArray(exposureGroupName, "colPix", colPix.data(),  colPix.size());
            hdf5File.writeArray(exposureGroupName, "flux",   flux.data(),    flux.size());            
        }
    }

    // Write the total sky background flux values [photons/pixel/exposure] to HDF5 in a custom group

    hdf5File.writeArray("Background/", "skyBackground", skyBackgroundValues.data(), skyBackgroundValues.size());

}












/**
 * \brief Configure the Camera object using the ConfigurationParameters
 * 
 * \param configParam: the configuration parameters 
 */

void Camera::configure(ConfigurationParameters &configParam)
{
    plateScale             = configParam.getDouble("Camera/PlateScale");
    focalLength            = configParam.getDouble("Camera/FocalLength") * 1000;                  // [m] -> [mm]
    throughputBandwidth    = configParam.getDouble("Camera/ThroughputBandwidth");
    throughputLambdaC      = configParam.getDouble("Camera/ThroughputLambdaC");

    polynomialType         = configParam.getString("Camera/FieldDistortion/Type");
    polynomialDegree       = configParam.getInteger("Camera/FieldDistortion/Degree");
    polynomialCoefficients = configParam.getDoubleVector("Camera/FieldDistortion/Coefficients");
    inversePolynomialCoefficients = configParam.getDoubleVector("Camera/FieldDistortion/InverseCoefficients");

    includeFieldDistortion = configParam.getBoolean("Camera/IncludeFieldDistortion");

    fluxOfV0Star           = configParam.getDouble("ObservingParameters/Fluxm0");                 // [phot/s/m^2/nm]
    userGivenSkyBackground = configParam.getDouble("ObservingParameters/SkyBackground");          // [phot/pix/s]
}










/** 
 * \brief      Specify the type of fit function used to fit the distortion
 *
 * \param      polynomial         The polynomial that describes the distortion
 * \param      inversePolynomial  The inverse polynomial
 */

void Camera::setDistortionPolynomial(Polynomial1D &polynomial, Polynomial1D &inversePolynomial)
{
    this->polynomial = polynomial;
    this->inversePolynomial = inversePolynomial;
}













/**
 * \brief      Expose the subField to the Sky, i.e. add flux to the detector,
 *             add the sky background, and convolve with the PSF.
 *
 * \param[in]  detector      the Detector class
 * \param[in]  startTime     start time of the exposure [seconds]
 * \param[in]  exposureTime  duration of one exposure [seconds]
 */

void Camera::exposeDetector(Detector &detector, double startTime, double exposureTime)
{
    // Get the focal plane coordinates of the center and the corners of the subfield (in [mm]).
    // To compute the diagonal length of the subfield, we only need the lower left (X00, Y00)
    // and the upper right (X11, Y11) corner of the subfield.

    double centerXmm, centerYmm;
    tie(centerXmm, centerYmm) = detector.getFocalPlaneCoordinatesOfSubfieldCenter();

    double corner00Xmm, corner00Ymm, corner11Xmm, corner11Ymm, dummy;
    tie(corner00Xmm, corner00Ymm, dummy, dummy, corner11Xmm, corner11Ymm, dummy, dummy) = detector.getFocalPlaneCoordinatesOfSubfieldCorners();

    // Convert the undistorted [mm] to distorted [mm] focal plane coordinates

    if (includeFieldDistortion)
    {
        Log.info("Camera: including field distortion");

        tie(centerXmm, centerYmm) = distortedToUndistortedFocalPlaneCoordinates(centerXmm, centerYmm);
        tie(corner00Xmm, corner00Ymm) = distortedToUndistortedFocalPlaneCoordinates(corner00Xmm, corner00Ymm);
        tie(corner11Xmm, corner11Ymm) = distortedToUndistortedFocalPlaneCoordinates(corner11Xmm, corner11Ymm);
    }

    Log.debug("Camera: center of subfield at (Xmm, Ymm) = (" + to_string(centerXmm) + ", " + to_string(centerYmm) + ") mm");
    Log.debug("Camera: lower left corner of subfield at (Xmm, Ymm) = (" + to_string(corner00Xmm) + ", " + to_string(corner00Ymm) + ") mm");
    Log.debug("Camera: upper right corner of subfield at (Xmm, Ymm) = (" + to_string(corner11Xmm) + ", " + to_string(corner11Ymm) + ") mm");


    // Convert the focal plane coordinates [mm] to (alpha, delta) equatorial sky coordinates [rad]

    double centerRA, centerDec;
    tie(centerRA, centerDec) = newFocalPlaneToSkyCoordinates(centerXmm, centerYmm);

    double corner00RA, corner00Dec;
    tie(corner00RA, corner00Dec) = newFocalPlaneToSkyCoordinates(corner00Xmm, corner00Ymm);

    double corner11RA, corner11Dec;
    tie(corner11RA, corner11Dec) = newFocalPlaneToSkyCoordinates(corner11Xmm, corner11Ymm);

    Log.debug("Camera: center of subfield at (alpha, delta) = (" + to_string(rad2deg(centerRA)) + ", " + to_string(rad2deg(centerDec)) + ") deg");
    Log.debug("Camera: lower left corner of subfield at (alpha, delta) = (" + to_string(rad2deg(corner00RA)) + ", " + to_string(rad2deg(corner00Dec)) + ") deg");
    Log.debug("Camera: upper right corner of subfield at (alpha, delta) = (" + to_string(rad2deg(corner11RA)) + ", " + to_string(rad2deg(corner11Dec)) + ") deg");


    // Compute the angular distance on the sky between the lower left and the upper right corner
    // of the subfield, to estimate the "radius" of the subfield.

    SkyCoordinates skyCoordinates00(corner00RA, corner00Dec, Angle::radians);
    SkyCoordinates skyCoordinates11(corner11RA, corner11Dec, Angle::radians);
    
    double radius = angularDistanceBetween(skyCoordinates00, skyCoordinates11, Angle::radians) / 2.0;

    Log.debug("Camera: semi-diagonal of subfield = " + to_string(rad2deg(radius)) + " deg");


    // Get a catalog of stars that fall on the subfield. Take the radius a bit larger so that the 
    // queried area includes possible small shifts of the projected subfield because of jitter.

    auto starCatalog = sky.getStarsWithinRadiusFrom(centerRA, centerDec, radius * 1.1, Angle::radians);

    Log.info("Camera: Found " + to_string(starCatalog.size()) + " stars on and near the subfield");  


    // If the telescope and/or platform show small variations (e.g. due to jitter) during the exposure,
    // the exposure time is split up in many small intervals, to track the effect of these variations
    // on the exposure. The largest time interval for which the variations can still be reliably tracked
    // is called the heartbeatInterval. The time step used should be either the heartbeat interval or the
    // expsosure time whatever is smallest. 

    double timeStep = min(telescope.getHeartbeatInterval(), exposureTime);


    // Later on we will have to convert from magnitudes to fluxes. Precompute a constant prefactor.
    // fluxOfV0Star is the photon flux [photons/s/m^2/nm] for a V=0 G2V-star.
    // Units of fluxFactor: [photons/s]
  
    const double fluxFactor = fluxOfV0Star * throughputBandwidth * telescope.getTransmissionEfficiency() * telescope.getLightCollectingArea(); 

    // Update the internal clock

    internalTime = startTime;

    // Take the flux of point sources (stars) into account.
    // Break up the exposure time in small intervals (hearbeat intervals) to track jitter while exposing.

    while (internalTime < startTime + exposureTime)
    {
        // Let the telescope pointing evolve over a small time interval

        telescope.updatePointingCoordinates(internalTime);

        // Loop over all stars in the catalog, and add their flux to the subfield

        for (int n = 0; n < starCatalog.size(); n++)
        {
            // Get the focal plane coordinates (in [mm]) of this particular star
            
            auto star = starCatalog[n];
            double Xmm, Ymm;
            tie(Xmm, Ymm) = newSkyToFocalPlaneCoordinates(star.RA, star.dec);

            // If required, include field distortion

            if (includeFieldDistortion)
            {
                tie(Xmm, Ymm) = undistortedToDistortedFocalPlaneCoordinates(Xmm, Ymm);
            }

            // Compute the flux [photons] of this star
            // Photons are always an integer number, so round down.

            double flux = floor(fluxFactor * pow(10.0, -0.4 * star.Vmag) * timeStep);

            // Let the detector add the flux to the appropriate pixel. 
            // Detector.flux() returns the pixel coordinates to which the flux was added.

            bool isInSubfield;
            double rowPix, colPix;    // subfield (not CCD) pixel coordinates

            tie(isInSubfield, rowPix, colPix) = detector.addFlux(Xmm, Ymm, flux);

            // If the star is indeed in the subfield, collect the following information to later write to HDF5
            //    1) average (Xmm, Ymm) coordinates of the star during the exposure                   [mm]
            //    2) average (row, col) pixel coordinates of the star on the CCD during the exposure  [pix]
            //    3) the total number of photons gathered of this star during the exposure            [photons]
            //    4) the total number of times that the star was in the subfield during the exposure
            //
            // Note: Due to jitter, the star can move in and out the subfield during the exposure

            if (isInSubfield)
            {
                // If this is the first time we encounter this startTime, initialise the information

                if (detectedStarInfo.find(startTime) == detectedStarInfo.end())
                {
                    detectedStarInfo[startTime][star.ID] = {{Xmm, Ymm, rowPix, colPix, flux, 1.0}};
                }
                else
                {
                    // If this is the first time that we encounter this star ID associated with this startTime,
                    // initialise the information. If not, just update the info.

                    if (detectedStarInfo[startTime].find(star.ID) == detectedStarInfo[startTime].end())
                    {
                        detectedStarInfo[startTime][star.ID] = {{Xmm, Ymm, rowPix, colPix, flux, 1.0}};
                    }
                    else
                    {
                        detectedStarInfo[startTime][star.ID][0] += Xmm;      // Will be used to compute average Xmm during the exposure
                        detectedStarInfo[startTime][star.ID][1] += Ymm;      // Will be used to compute average Ymm during the exposure
                        detectedStarInfo[startTime][star.ID][2] += rowPix;   // Will be used to compute average pixel row during the exposure
                        detectedStarInfo[startTime][star.ID][3] += colPix;   // Will be used to compute average pixel column during the exposure
                        detectedStarInfo[startTime][star.ID][4] += flux;     // Total flux
                        detectedStarInfo[startTime][star.ID][5] += 1;        // # of times a star was on the subfield during an exposure 
                    }
                }
            }
        }

        Log.debug("Camera: Incremented flux of stars in subfield");


        // Update the clock. Normally with 'timeStep', but if adding timeStep would overstep
        // the total exposure time, take the small rest time instead.

        timeStep = min(timeStep, startTime + exposureTime - internalTime);
        internalTime += timeStep;

    }

    // Take the flux of the stellar background and the zodiacal light into account. Use one value for the entire subfield.
    // A negative value for the user given sky background value [phot/pix/s] signals that we should compute it ourselves.
    // Note: - The output of sky.zodiacalFlux() is in [J s^{-1} m^{-2} sr^{-1} m^{-1}]
    //       - As wavelength range we take the entire throughput band.
    //       - Photons are always an integer number, thus round down.

    totalSkyBackground = 0.0;


    if (userGivenSkyBackground < 0.0)
    {
        const double energyOfOnePhoton = Constants::CLIGHT * Constants::HPLANCK / (throughputLambdaC * 1.e-9);                // [J]
        const double lambda1 = (throughputLambdaC - throughputBandwidth/2.0) * 1.e-9;                                         // [m]
        const double lambda2 = (throughputLambdaC + throughputBandwidth/2.0) * 1.e-9;                                         // [m]
    
        const double zodiacalFlux = sky.zodiacalFlux(centerRA, centerDec, lambda1, lambda2)                                                      // [phot/exposure]
                                    * exposureTime * telescope.getTransmissionEfficiency() * telescope.getLightCollectingArea()
                                    * detector.getSolidAngleOfOnePixel(plateScale) / energyOfOnePhoton; 

        const double stellarBackgroundFlux = sky.stellarBackgroundFlux(centerRA, centerDec, lambda1, lambda2)                                    // [phot/exposure]
                                             * exposureTime * telescope.getTransmissionEfficiency() * telescope.getLightCollectingArea()
                                             * detector.getSolidAngleOfOnePixel(plateScale) / energyOfOnePhoton;      


        totalSkyBackground = floor(zodiacalFlux + stellarBackgroundFlux);
        detector.addFlux(totalSkyBackground);

        Log.debug("Camera: zodiacal flux level in subfield = " + to_string(zodiacalFlux) + " photons/pixel/exposure");
        Log.debug("Camera: stellar background flux level in subfield = " + to_string(stellarBackgroundFlux) + " photons/pixel/exposure");
    }
    else
    {
        totalSkyBackground = floor(userGivenSkyBackground * exposureTime);
        detector.addFlux(totalSkyBackground);

        Log.debug("Camera: user-given sky background flux = " + to_string(userGivenSkyBackground * exposureTime) + " photons/pixel/exposure");
    }

    // Save the sky background value that we added. [photons/pix/exposure]

    skyBackgroundValues.push_back(totalSkyBackground);

    // Convolve with the point spread function
    // Detector was given the proper PSF in Simulation::run().

    detector.convolveWithPsf();

    return;
}












/**
 * \brief      Select the PSF for the given focal plane coordinates
 *
 * \details    This method selects, rotates and rebins the PSF.
 *
 * \param[in]  xFPmm             x-coordinate in the FP' reference frame [mm]
 * \param[in]  yFPmm             y-coordinate in the FP' reference frame [mm]
 * \param[in]  targetSubPixels   the number of subpixels per pixels in the detector subfield
 * \param[in]  orientationAngle  the orientation of the CCD wrt focal plane orientation (counter clockwise) [rad]
 *
 * \return     the psfMap that was selected, rotated and rebinned
 */
arma::fmat Camera::getRebinnedPsfForFocalPlaneCoordinates(double xFPmm, double yFPmm, unsigned int targetSubPixels, double orientationAngle)
{
    arma::fmat psfMap;

    // Get the 'user specified' angular distance to the optical axis from the psf.
    // If the user didn't specify an angular distance, calculate it from the given
    // focal plane coordinates.

    double radius = psf->getRequestedDistanceToOpticalAxis();

    if (radius < 0.0)
    {
        radius = getGnomonicRadialDistanceFromOpticalAxis(xFPmm, yFPmm);
    }

    psf->select(radius);


    // Get the 'user specified' orientation angle from the psf.
    // if the user didn't specify a rotation angle, calculate it 
    // from the given focal plane coordinates.

    double angle = psf->getRequestedRotationAngle();

    if (angle < 0.0)
    {
        angle = atan2(yFPmm, xFPmm);
    }

    //  Compensate for the orientation of the CCD wrt focal plane orientation.

    angle -= orientationAngle;

    psf->rotate(angle);


    // Rebin the psfMap to the number of sub-pixels per pixel used for the Detector

    psfMap = psf->rebinToSubPixels(targetSubPixels);

    return psfMap;
}


















/**
 * @brief      Calculate the gnomonic radial distance with respect to the optical axis in the focal plane
 *
 * @param[in]  xFPprime  Focal plane x-coordinate [mm]
 * @param[in]  yFPprime  Focal plane y-coordinate [mm]
 *
 * @return     the angular distance of the star w.r.t. the optical axis [rad]
 * 
 */
double Camera::getGnomonicRadialDistanceFromOpticalAxis(double xFPprime, double yFPprime)
{
    const double tanx = xFPprime / focalLength;
    const double tany = yFPprime / focalLength;

    double angularDistance = acos(1.0/sqrt(1.0 + tanx*tanx + tany*tany));

    // Take care that the angle is between [0, 2*PI]

    if (angularDistance < 0.0)
    {
        angularDistance += 2.0 * Constants::PI;
    }

    if (angularDistance > 2.0 * Constants::PI)
    {
        angularDistance -= 2.0 * Constants::PI;
    }

    // That's it!

    return angularDistance;
}












/**
 * \brief      Computes the (x,y) coordinates in the focal plane, of a star with given equatorial 
 *             sky coordinates, assuming a pinhole camera.
 *             
 * \details    The transformation is with respect to the current pointing coordinates of the telescope.
 *
 * \param raStar       Right ascension of the star [rad]
 * \param decStar      Declination of the star [rad]
 *
 * return pair (x,y):  Cartesian coordinate of the projected star in the focal plane in the FP-prime system [mm]
 */

pair<double, double> Camera::skyToFocalPlaneCoordinates(double raStar, double decStar)
{
    // Get the current equatorial coordinates of the optical axis [rad]

    double raOpticalAxis, decOpticalAxis;
    tie(raOpticalAxis, decOpticalAxis) = telescope.getCurrentPointingCoordinates();

    // Get the current focal plane orientation

    const double focalPlaneAngle = telescope.getCurrentFocalPlaneOrientation();

    return skyToFocalPlaneCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle);
}














/**
 * \brief      Computes the (x,y) coordinates in the focal plane, of a star with given equatorial 
 *             sky coordinates, assuming a pinhole camera.
 *             
 * \details    The transformation is with respect to the given pointing coordinates of the telescope.
 *
 * \param raStar          Right ascension of the star [rad]
 * \param decStar         Declination of the star [rad]
 * \param raOpticalAxis   Right ascension of the optical axis [rad]
 * \param decOpticalAxis  Declination of the optical axis [rad]
 * \param focalPlaneAngle The Focal Plane Orientation [rad]
 *
 * return pair (x,y):  Cartesian coordinate of the projected star in the focal plane in the FP-prime system [mm]
 */

pair<double, double> Camera::skyToFocalPlaneCoordinates(double raStar, double decStar, double raOpticalAxis, double decOpticalAxis, double focalPlaneAngle)
{

    // Convert the equatorial sky coordinate of the star to equatorial cartesian coordinates on the unit sphere

    arma::vec vecEQ = {cos(decStar) * cos(raStar), cos(decStar) * sin(raStar), sin(decStar)};

    // Convert the equatorial cartesian coordinates to focal plane cartesian coordinates.

    arma::mat rotMatrix1;
    rotMatrix1 <<  cos(raOpticalAxis) << sin(raOpticalAxis) <<  0 << arma::endr
               << -sin(raOpticalAxis) << cos(raOpticalAxis) <<  0 << arma::endr
               <<          0          <<          0         <<  1 << arma::endr;

    arma::mat rotMatrix2;
    rotMatrix2 << sin(decOpticalAxis) << 0 << -cos(decOpticalAxis) << arma::endr
               <<          0          << 1 <<          0           << arma::endr
               << cos(decOpticalAxis) << 0 <<  sin(decOpticalAxis) << arma::endr;

    arma::vec vecFP = rotMatrix2 * rotMatrix1 * vecEQ;

    // Take into account the projection effect of the pinhole camera
    // Note that the pinhole reverses the image, hence the minus signs.
    // The focalLength is assumed to be in [mm] so that the (xFP, yFP) coordinates are also in [mm].

    const double xFP = -focalLength * vecFP(0)/vecFP(2);
    const double yFP = -focalLength * vecFP(1)/vecFP(2);

    // Convert the FP coordinates into FP' coordinates 

    const double xFPprime =  xFP * cos(focalPlaneAngle) + yFP * sin(focalPlaneAngle);
    const double yFPprime = -xFP * sin(focalPlaneAngle) + yFP * cos(focalPlaneAngle);

    // That's it!

    return make_pair(xFPprime, yFPprime);
}















/**
 * \brief Compute the equatorial sky coordinates of a star which has the given focal plane (FP') coordinates (x,y),
 *        assuming a pinhole camera model
 *        
 * \param xFPprime     Focal plane x-coordinate in the FP-prime system [mm]
 * \param yFPprime     Focal plane y-coordinate in the FP-prime system [mm]
 *
 * \return (alpha, delta)  Equatorial coordinates (RA & Dec) of the star [rad]
 */

pair<double, double> Camera::focalPlaneToSkyCoordinates(double xFPprime, double yFPprime)
{    
    // Get the current equatorial coordinates of the optical axis [rad]

    double raOpticalAxis, decOpticalAxis;
    tie(raOpticalAxis, decOpticalAxis) = telescope.getCurrentPointingCoordinates();

    // Get the current focal plane orientation

    const double focalPlaneAngle = telescope.getCurrentFocalPlaneOrientation();

    // Convert the FP' coordinates in FP coordinates

    const double xFP =  xFPprime * cos(focalPlaneAngle) - yFPprime * sin(focalPlaneAngle);
    const double yFP =  xFPprime * sin(focalPlaneAngle) + yFPprime * cos(focalPlaneAngle);

    // Undo the reverse-image projection effect of the pinhole.
    // Both the focalLength and the (xFP, yFP) coordinates are assumed to be in [mm].

    arma::vec vecFP = {-xFP/focalLength, -yFP/focalLength, 1.0};

    // Convert the focal plane cartesian coordinates to equatorial cartesian coordinates. 

    arma::mat rotMatrix1;
    rotMatrix1 <<  sin(decOpticalAxis) << 0 << cos(decOpticalAxis) << arma::endr
               <<          0           << 1 <<          0          << arma::endr
               << -cos(decOpticalAxis) << 0 << sin(decOpticalAxis) << arma::endr;

    arma::mat rotMatrix2;
    rotMatrix2 << cos(raOpticalAxis) << -sin(raOpticalAxis) <<  0 << arma::endr
               << sin(raOpticalAxis) <<  cos(raOpticalAxis) <<  0 << arma::endr
               <<          0         <<          0          <<  1 << arma::endr;

    arma::vec vecEQ = rotMatrix2 * rotMatrix1 * vecFP;


    // Convert the cartesian equatorial coordinates to equatorial sky coordinates

    const double norm = sqrt(vecEQ(0)*vecEQ(0) + vecEQ(1)*vecEQ(1) + vecEQ(2)*vecEQ(2)); 
    double decStar = Constants::PI/2.0 - acos(vecEQ(2)/norm);
    double raStar = atan2(vecEQ(1), vecEQ(0));

    if (raStar < 0.0)
    {
        raStar += 2.*Constants::PI;
    }

    // That's it!

    return make_pair(raStar, decStar);
}










/**
 * \brief      Convert polar coordinates to cartesian coordinates
 *
 * \param[in]  distance  distance from the pole (reference point) [mm]
 * \param[in]  angle     angle counter-clockwise from the x-axis [rad]
 *
 * \return     (xFPmm, yFPmm) Cartesian coordinates in the focal plane [mm]
 */
pair<double, double> Camera::polarToCartesianFocalPlaneCoordinates(double distance, double angle)
{
    double xFPmm = cos(angle) * distance;
    double yFPmm = sin(angle) * distance;

    return make_pair(xFPmm, yFPmm);
}





/**
 * \brief      Convert cartesian coordinates to polar coordinates
 *
 * \param[in]  xFPmm  x-axis cartesian coordinate in the focal plane [mm]
 * \param[in]  yFPmm  y-axis cartesian coordinate in the focal plane [mm]
 *
 * \return     (distance, angle) polar coordinates in the focal plane
 */
pair<double, double> Camera::cartesianToPolarFocalPlaneCoordinates(double xFPmm, double yFPmm)
{
    double angle = atan2(yFPmm, xFPmm);
    double distance = sqrt(xFPmm * xFPmm + yFPmm * yFPmm);

    return make_pair(distance, angle);
}












/**
 * @brief      Convert from undistorted to distorted focal plane coordinates
 *
 * @param[in]  xFPmm  Undistorted focal plane x-coordinate [mm]
 * @param[in]  yFPmm  Undistorted focal plane y-coordinate [mm]
 *
 * @return     (xFPdist, yFPdist) distorted x and y coordinates [mm]
 */
pair<double, double> Camera::undistortedToDistortedFocalPlaneCoordinates(double xFPmm, double yFPmm)
{
    double alpha = atan2(yFPmm, xFPmm);  // [radians]
    
    double rFP = sqrt(xFPmm * xFPmm + yFPmm * yFPmm);
    double rFPdist = polynomial(rFP);
    
    double xFPdist = cos(alpha) * rFPdist;
    double yFPdist = sin(alpha) * rFPdist;
    
    return make_pair(xFPdist, yFPdist);
}







/**
 * @brief      Convert from distorted to undistorted focal plane coordinates
 *
 * @param[in]  xFPdist  Distorted focal plane x-coordinate [mm]
 * @param[in]  yFPdist  DIstorted focal plane y-coordinate [mm]
 *
 * @return     (xFPmm, yFPmm) distorted x and y coordinates [mm]
 */
pair<double, double> Camera::distortedToUndistortedFocalPlaneCoordinates(double xFPdist, double yFPdist)
{
    double alpha = atan2(yFPdist, xFPdist);  // [radians]
    
    double rFP   = sqrt(xFPdist * xFPdist + yFPdist * yFPdist);
    double rFPmm = inversePolynomial(rFP);
    
    double xFPmm = cos(alpha) * rFPmm;
    double yFPmm = sin(alpha) * rFPmm;
    
    return make_pair(xFPmm, yFPmm);
}






/**
 * @brief   Returns the total sky background, expressed in photons / pixel / exposure.
 *
 * @ return  Total sky background [photons / pixel / exposure].
 */
double Camera::getTotalSkyBackground()
{
	return totalSkyBackground;
}
