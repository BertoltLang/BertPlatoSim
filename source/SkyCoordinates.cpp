
#include "SkyCoordinates.h"


/**
 * \brief constructor
 * 
 * \param RA: right ascension 
 * \param decl: declination
 * \param angleUnit: Angle::radians if input angles are in radians, Angle::degrees if in degrees
 *                   
 */

SkyCoordinates::SkyCoordinates(double RA, double decl, Unit angleUnit)
: RA(RA), decl(decl), obliquity(0.409087723), inclGalPlane(1.0925761117484503), 
  alphaN(4.926191813753995), l0(0.5759586531581287)
{
    this->RA = RA / angleUnit;
    this->decl = decl / angleUnit;
}






/**
 * \brief destructor
 * 
 */

SkyCoordinates::~SkyCoordinates()
{
} 








/**
 * \brief Return the galactic coordinates l and b
 * 
 * \details The conversion is done for epoch B1950. 
 *          Example: double l,b;
 *                   tie(l,b) = point.toGalactic(Angle::degrees);
 * 
 * \param angleUnit: Angle::radians if output angles should be in radians, Angle::degrees if in degrees
 * \return (l, b)
 */

pair<double, double> SkyCoordinates::toGalactic(Unit angleUnit)
{
    double b,l;
    equatorial2galactic(RA, decl, l, b);
    return make_pair(l * angleUnit, b * angleUnit);
}








/**
 * \brief return the ecliptic coordinates lambda and beta
 * 
 * \details Example: double l,b;
 *                   tie(l,b) = point.toGalactic(Angle::degrees);
 * 
 * \param angleUnit: Angle::radians if output angles should be in radians, Angle::degrees if in degrees
 * \return (lambda, beta)
 */

pair<double, double> SkyCoordinates::toEcliptic(Unit angleUnit)
{
    double beta,lambda;
    equatorial2ecliptic(RA, decl, lambda, beta);
    return make_pair(lambda * angleUnit, beta * angleUnit);
}







/**
 * \brief Convert equatorial coordinates (alpha, delta) into ecliptic coordinates (lambda, beta)
 * 
 * \note  Epoch=2000.0
 * 
 * \param alpha[in]:   equtorial right ascension [rad]
 * \param delta[in]:   equatorial declination [rad]
 * \param lambda[out]: ecliptic longitude [rad]
 * \param beta[out]:   ecliptic latitude [rad]
 */

void SkyCoordinates::equatorial2ecliptic(const double alpha, const double delta, double &lambda, double &beta)
{
    double sinbeta = sin(delta) * cos(obliquity) - cos(delta) * sin(obliquity) * sin(alpha);
    beta = asin(sinbeta);
    double cosbeta = cos(beta);

    if (cosbeta == 0.0)
    {
        Log.error("equatorial2ecliptic: pointing to ecliptic pole.");
        exit (1);
    }

    double sinlambda = (sin(delta) * sin(obliquity) + cos(delta) * cos(obliquity) * sin(alpha)) / cosbeta;
    double coslambda = cos(alpha) * cos(delta) / cosbeta;

    lambda = atan2(sinlambda, coslambda);

    if (lambda < 0.0) lambda += 2.0 * 3.141592653589793;

}








/**
 * \brief Convert equatorial coordinates (alpha, delta) into galactic coordinates (l, b)
 * 
 * \details The conversion is done for epoch B1950. 
 * 
 * \param alpha[in]: equtorial right ascension [rad]
 * \param delta[in]: equatorial declination [rad]
 * \param l[out]:    galactic longitude [rad]
 * \param b[out]:    galactic latitude [rad]
 */

void SkyCoordinates::equatorial2galactic(const double alpha, const double delta, double &l, double &b)
{
    double sinb = sin(delta) * cos(inclGalPlane) - cos(delta) * sin(alpha - alphaN) * sin(inclGalPlane);
    b = asin(sinb);
    double cosb = cos(b);

    if (cosb == 0.0)
    {
        Log.error("equatorial2galactic: pointing to galactic pole.");
        exit (1);
    }

    double sine = (sin(delta) * sin(inclGalPlane) + cos(delta) * cos(inclGalPlane) * sin(alpha-alphaN)) / cosb;
    double cosine = cos(alpha-alphaN) * cos(delta) / cosb;

    l = atan2(sine, cosine) + l0;

    if (l < 0.0) l += 2.0 * 3.141592653589793;

} 










/**
 * \brief Compute the angular separation on a great-circle between two coordinate pairs
 * 
 * \details The haversine formula is used to get better numerical accuracy.
 * 
 * \param coordinates1: coordinates of point1
 * \param coordinates2: coordinates of point2
 * \param angleUnit: Angle::radians if output angle should be in radians, Angle::degrees if in degrees
 *   
 * \return angular separation [rad]
 */

double angularDistanceBetween(SkyCoordinates &coordinates1, SkyCoordinates &coordinates2, Unit outputAngleUnit)
{
    const double sinHalfDeltaLong = sin((coordinates1.RA - coordinates2.RA)/2.0);
    const double sinHalfDeltaLat = sin((coordinates1.decl - coordinates2.decl)/2.0);
    const double a = sinHalfDeltaLat*sinHalfDeltaLat + cos(coordinates1.decl) * cos(coordinates2.decl) * sinHalfDeltaLong*sinHalfDeltaLong;

    const double angle = 2.0 * atan2(sqrt(a), sqrt(1.0-a)) * outputAngleUnit;

    return angle;
}











/**
 * \brief Compute the angular distance of for each of the stars with coordinates (RA, dec) to a reference point (RA0, dec0).
 *        
 * \details The haversine formula is used to get better numerical accuracy.
 * 
 * \param RA0              Right ascension of the reference point [rad]
 * \param dec0             Declination of the reference point [rad]
 * \param RA               Right ascension of each of the stars [rad]
 * \param dec              Declination of each of the stars [rad]
 * \param outputAngleUnit  anglesAngle::radians if output angle should be in radians, Angle::degrees if in degrees
 *   
 * \return           Angular separation [unit: see angleUnit]
 */

vector<double> angularDistanceBetween(const double RA0, const double dec0, const vector<double> &RA, const vector<double> &dec, Unit outputAngleUnit)
{
    vector<double> angularDistance(RA.size());

    for (long n = 0; n < RA.size(); ++n)
    {
        const double sinHalfDeltaLong = sin((RA0 - RA[n])/2.0);
        const double sinHalfDeltaLat = sin((dec0 - dec[n])/2.0);
        const double a = sinHalfDeltaLat*sinHalfDeltaLat + cos(dec0) * cos(dec[n]) * sinHalfDeltaLong*sinHalfDeltaLong;

        angularDistance[n] = 2.0 * atan2(sqrt(a), sqrt(1.0-a)) * outputAngleUnit;
    }

    return angularDistance;
}
