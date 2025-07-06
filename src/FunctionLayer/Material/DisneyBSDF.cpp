#include "DisneyBSDF.h"
#include "BxDF/MicrofacetDistribution.h"
// #include "BxDF/Fresnel.h"
// #include  "FunctionLayer/Texture/TextureFactory.h"
// #include "FunctionLayer/Distribution/Distribution.h"
// #include "FastMath.h"

#include <variant>

//for debug purpose
static bool onlyDiffuse = false;
static bool onlyMetal = false;
static bool onlyGlass = false;
static bool onlySheen = false;
static bool onlyClearCoat = false;

static Vector2f  getAnisAlpha(double  roughness, double anis) {
    double aspect = sqrt(1 - 0.9 * anis);
    float alphaX = std::max(0.0001, pow(roughness, 2) / aspect);
    float alphaY = std::max(0.0001, pow(roughness, 2) * aspect);
    return { alphaX,alphaY };
}


static  double smithMaskingGTR2(const Vector3f& v_local, double roughness) {
    double alpha = roughness * roughness;
    double a2 = alpha * alpha;
    Vector3f v2 = v_local * v_local;
    double Lambda = (-1 + sqrt(1 + (v2.x() * a2 + v2.y() * a2) / v2.z())) / 2;
    return 1 / (1 + Lambda);
}

struct DisneyDiffuse {
    const Spectrum  baseColor;
    const double roughness;
    const double subsurface;
};
struct DisneyMetal {
    const Spectrum  baseColor;
    Vector2f alpha;
    GGXDistribution distrib;
    DisneyMetal(const Spectrum& baseColor, double roughness, double anis) :baseColor(baseColor), distrib(true) {
        roughness = std::clamp(roughness, double(0.01), double(1));
        alpha = getAnisAlpha(roughness, anis);
    }
};
struct DisneyClearCoat {
    double clearCoatGloss;
};
struct DisneyGlass {
    const Spectrum  baseColor;
    Vector2f alpha;
    double eta;
    GGXDistribution distrib;
    DisneyGlass(const Spectrum& baseColor, double roughness, double anis, double eta) : baseColor(baseColor), eta(eta),
        distrib(true) {
        roughness = std::clamp(roughness, double(0.01), double(1));
        alpha = getAnisAlpha(roughness, anis);
    }
};
struct DisneySheen {
    const Spectrum  baseColor;
    double sheenTint;
};

using DisneyBXDF = std::variant<DisneyDiffuse,
    DisneyMetal,
    DisneyClearCoat,
    DisneyGlass,
    DisneySheen>;


struct evalDisneyBXDFOP {
    Spectrum operator()(const DisneyDiffuse& disneyBXDF);
    Spectrum operator()(const DisneyMetal& disneyBXDF);
    Spectrum operator()(const DisneyClearCoat& disneyBXDF);
    Spectrum operator()(const DisneyGlass& disneyBXDF);
    Spectrum operator()(const DisneySheen& disneyBXDF);
    const Vector3f& out;
    const Vector3f& in;
};

struct pdfDisneyBXDFOP {
    double operator()(const DisneyDiffuse& disneyBXDF);
    double operator()(const DisneyMetal& disneyBXDF);
    double operator()(const DisneyClearCoat& disneyBXDF);
    double operator()(const DisneyGlass& disneyBXDF);
    double operator()(const DisneySheen& disneyBXDF);
    const Vector3f& out;
    const Vector3f& in;
};

struct sampleDisneyBXDFOP {
    BSDFSampleResult operator()(const DisneyDiffuse& disneyBXDF);
    BSDFSampleResult operator()(const DisneyMetal& disneyBXDF);
    BSDFSampleResult operator()(const DisneyClearCoat& disneyBXDF);
    BSDFSampleResult operator()(const DisneyGlass& disneyBXDF);
    BSDFSampleResult operator()(const DisneySheen& disneyBXDF);
    const Vector3f& out;
    const Vector2f& sample;
};

///-----------------------------Disney Diffuse Functions Begin -----------------------------///

Spectrum evalDisneyBXDFOP::operator ()(const DisneyDiffuse& disneyBXDF) {
    if (in.z() < 0 || out.z() < 0)
        return Spectrum(0);
    Vector3f wh = normalize(in + out);

    const double roughness = disneyBXDF.roughness;
    const double  subSurfaceFactor = disneyBXDF.subsurface;
    const Spectrum  baseColor = disneyBXDF.baseColor / M_PI;

    double FD90 = 0.5 + 2 * roughness * abs(dot(wh, in)) * abs(dot(wh, in));
    auto FD = [&](const Vector3f& w) {
        return 1 + (FD90 - 1) * pow((1 - AbsCosTheta(w)), 5);
        };
    Spectrum baseDiffuse = baseColor * FD(out) * FD(in);

    double FSS90 = roughness * abs(dot(wh, in));
    auto FSS = [&](const Vector3f& w) {
        return 1 + (FSS90 - 1) * pow((1 - AbsCosTheta(w)), 5);
        };
    double scaleFactor =
        FSS(in) * FSS(out) * (1 / (AbsCosTheta(in) + AbsCosTheta(out)) - 0.5) + 0.5;
    Spectrum subsurface = 1.25 * baseColor * scaleFactor;

    return ((1 - subSurfaceFactor) * baseDiffuse + subSurfaceFactor * subsurface);
}

double pdfDisneyBXDFOP::operator ()(const DisneyDiffuse& disneyBXDF) {
    if (in.z() < 0 || out.z() < 0)
        return 0;
    return std::max(in.z(), 0.0f) / M_PI;
}

BSDFSampleResult sampleDisneyBXDFOP::operator ()(const DisneyDiffuse& disneyBXDF) {
    BSDFSampleResult result;
    Vector3f in = SquareToCosineHemisphere(sample);
    result.wi = in;
    result.type = BXDFType(BXDF_REFLECTION | BXDF_GLOSSY);
    return result;
}

///-----------------------------Disney Diffuse Functions End -----------------------------///

///-----------------------------Disney Metal Functions Begin -----------------------------///

Spectrum evalDisneyBXDFOP::operator ()(const DisneyMetal& disneyBXDF) {
    if (in.z() < 0 || out.z() < 0)
        return Spectrum(0);
    Vector3f wh = normalize(in + out);
    Spectrum baseColor = disneyBXDF.baseColor;
    Spectrum F = baseColor + (Spectrum(1.0) - baseColor) * pow((1 - AbsCosTheta(wh)), 5);
    Vector2f alpha = disneyBXDF.alpha;
    double G = disneyBXDF.distrib.G1(in, alpha) * disneyBXDF.distrib.G1(out, alpha);
    double D = disneyBXDF.distrib.D(wh, alpha);
    return F * D * G / (4 * AbsCosTheta(in) * AbsCosTheta(out));
}

double pdfDisneyBXDFOP::operator ()(const DisneyMetal& disneyBXDF) {
    if (in.z() < 0 || out.z() < 0) {
        return 0;
    }
    Vector3f wh = normalize(in + out);
    Vector2f alpha = disneyBXDF.alpha;
    return disneyBXDF.distrib.Pdf(out, wh, alpha) / (4 * dot(out, wh));
}


BSDFSampleResult sampleDisneyBXDFOP::operator ()(const DisneyMetal& disneyBXDF) {
    BSDFSampleResult result;
    if (CosTheta(out) < 0) {
        return {};
    }
    Vector2f alpha = disneyBXDF.alpha;
    Vector3f wh = disneyBXDF.distrib.Sample_wh(out, sample, alpha);
    Vector3f reflected = normalize(-out + 2 * dot(out, wh) * wh);
    result.wi = reflected;
    result.type = BXDFType(BXDF_REFLECTION | BXDF_GLOSSY);
    return result;
}

///-----------------------------Disney Metal Functions End -----------------------------///

///-----------------------------Disney ClearCoat Functions Begin -----------------------------///

Spectrum evalDisneyBXDFOP::operator ()(const DisneyClearCoat& disneyBXDF) {
    if (in.z() < 0 || out.z() < 0)
        return { 0 };
    Vector3f wh = normalize(in + out);
    double eta = 1.5;
    double R0 = pow((eta - 1) / (eta + 1), 2);
    double f = R0 + (1 - R0) * pow((1 - abs(dot(wh, in))), 5);
    Spectrum F = Spectrum(f);
    double roughness = disneyBXDF.clearCoatGloss;
    double alphaG = (1 - roughness) * 0.1 + roughness * 0.001;
    double D = (alphaG * alphaG - 1) / (M_PI * 2 * log(alphaG) * (1 + (alphaG * alphaG - 1) * pow(wh.z(), 2)));
    double G = smithMaskingGTR2(in, 0.25) * smithMaskingGTR2(out, 0.25);

    return F * D * G / (4 * abs(in.z() * out.z()));
}

double pdfDisneyBXDFOP::operator ()(const DisneyClearCoat& disneyBXDF) {
    if (CosTheta(out) < 0 || CosTheta(in) < 0) {
        return 0;
    }
    Vector3f wh = normalize(in + out);
    double cosTheta2 = wh.z() * wh.z();
    double roughness = disneyBXDF.clearCoatGloss;
    double alphaG = (1 - roughness) * 0.1 + roughness * 0.001;
    double alphaG2 = alphaG * alphaG;

    double  D = (alphaG2 - 1) / log(alphaG2) / M_PI / (1 + (alphaG2 - 1) * cosTheta2);
    return D * AbsCosTheta(wh) / (4 * abs(dot(wh, in)));
}

BSDFSampleResult sampleDisneyBXDFOP::operator ()(const DisneyClearCoat& disneyBXDF) {
    if (CosTheta(out) < 0) {
        return {};
    }
    BSDFSampleResult result;
    double roughness = disneyBXDF.clearCoatGloss;
    double alphaG = (1 - roughness) * 0.1 + roughness * 0.001;
    double alphaG2 = alphaG * alphaG;
    double cosTheta = sqrt((1 - pow(alphaG2, 1 - sample.x)) / (1 - alphaG2));
    double sinTheta = sqrt(1 - cosTheta * cosTheta);
    double phi = sample.y * M_PI * 2;
    Vector3f wh = Vector3f(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);

    Vector3f in = normalize(-out + 2 * (dot(out, wh)) * wh);
    result.wi = in;
    result.type = BXDFType(BXDF_GLOSSY | BXDF_REFLECTION);
    return result;
}

///-----------------------------Disney ClearCoat Functions End -----------------------------///

///-----------------------------Disney Glass Functions End -----------------------------///

Spectrum evalDisneyBXDFOP::operator ()(const DisneyGlass& disneyBXDF) {

    bool reflect = CosTheta(out) * CosTheta(in) > 0;

    Spectrum baseColor = disneyBXDF.baseColor;

    Vector2f alpha = disneyBXDF.alpha;

    double  eta = CosTheta(out) > 0 ? 1 / disneyBXDF.eta : disneyBXDF.eta;
    Vector3f wh;
    if (reflect) {
        wh = std::copysign(1, out.z()) * normalize(in + out);
    }
    else {
        wh = -normalize(in + out * eta);
    }
    double F = Fresnel::dielectricReflectance(1 / disneyBXDF.eta, dot(out, wh));
    double D = disneyBXDF.distrib.D(wh, alpha);
    double G = disneyBXDF.distrib.G1(out, alpha) * disneyBXDF.distrib.G1(in, alpha);

    if (reflect) {
        return baseColor * F * D * G / (4 * abs(in.z() * out.z()));
    }
    else {
        double whDotIn = dot(wh, in);
        double whDotOut = dot(wh, out);
        double sqrtDeom = eta * whDotOut + whDotIn;

        return  sqrt(baseColor) * (1 - F) * D * G * fm::abs(
            whDotIn * whDotOut /
            (in.z() * out.z() * sqrtDeom * sqrtDeom));
    }
}

double pdfDisneyBXDFOP::operator ()(const DisneyGlass& disneyBXDF) {
    bool reflect = CosTheta(out) * CosTheta(in) > 0;
    double eta = out.z() > 0 ? 1 / disneyBXDF.eta : disneyBXDF.eta;
    Vector3f wh;
    if (reflect) {
        wh = normalize(in + out) * std::copysign(1, out.z());
    }
    else {
        wh = -normalize(in + out * eta);
    }
    double  F = Fresnel::dielectricReflectance(1 / disneyBXDF.eta, dot(out, wh));
    auto alpha = disneyBXDF.alpha;
    double pdf;
    if (reflect) {
        pdf = F * disneyBXDF.distrib.Pdf(out, wh, disneyBXDF.alpha) / (4 * abs(dot(out, wh)));
    }
    else {

        double sqrtDenom = dot(out, wh) * eta + dot(in, wh);
        double dWhDWi =
            fm::abs(dot(in, wh)) / (sqrtDenom * sqrtDenom);
        pdf = disneyBXDF.distrib.Pdf(out, wh, alpha) * (1 - F) * dWhDWi;
    }
    return pdf;
}


BSDFSampleResult sampleDisneyBXDFOP::operator ()(const DisneyGlass& disneyBXDF) {
    BSDFSampleResult result;

    Vector3f wh = disneyBXDF.distrib.Sample_wh(out, sample, disneyBXDF.alpha);
    double whDotOut = dot(wh, out);
    double cosThetaT;
    double  F = Fresnel::dielectricReflectance(1 / disneyBXDF.eta, whDotOut, cosThetaT);
    double r = rand() % (10000 + 1) / (double)(10000 + 1);
    bool reflect = r < F;
    if (reflect) {
        result.wi = -out + 2 * dot(out, wh) * wh;
        result.type = BXDFType(BXDF_REFLECTION | BXDF_GLOSSY);
    }
    else {
        double eta = whDotOut < 0 ? disneyBXDF.eta : 1 / disneyBXDF.eta;
        Vector3f in = (eta * whDotOut - (whDotOut > 0 ? 1 : -1) * cosThetaT) * wh - eta * out;
        result.wi = in;
        result.type = BXDFType(BXDF_TRANSMISSION | BXDF_GLOSSY);
    }
    return result;
}

///-----------------------------Disney Glass Functions Begin-----------------------------///

///-----------------------------Disney Sheen Functions Begin -----------------------------///

Spectrum evalDisneyBXDFOP::operator ()(const DisneySheen& disneyBXDF) {
    if (in.z() < 0 || out.z() < 0) {
        return Spectrum(0);
    }
    Vector3f wh = normalize(in + out);
    double sheenTint = disneyBXDF.sheenTint;
    Spectrum baseColor = disneyBXDF.baseColor;
    Spectrum cTint = baseColor.luminance() > 0 ? baseColor / baseColor.luminance() : Spectrum(1.0);
    Spectrum cSheen = Spectrum(1 - sheenTint) + sheenTint * cTint;
    return cSheen * (pow(1 - abs(dot(wh, in)), 5));
}

double pdfDisneyBXDFOP::operator ()(const DisneySheen& disneyBXDF) {
    if (in.z() < 0 || out.z() < 0)
        return 0;
    return std::max(in.z(), 0.0f) / M_PI;
}

BSDFSampleResult sampleDisneyBXDFOP::operator ()(const DisneySheen& disneyBXDF) {
    BSDFSampleResult result;
    Vector3f in = SquareToCosineHemisphere(sample);
    result.wi = in;
    result.type = BSDFType(BXDF_REFLECTION | BXDF_GLOSSY);
    return result;
}

///-----------------------------Disney Sheen Functions End -----------------------------///


Spectrum evalDisneyBXDF(const DisneyBXDF& disneyBXDF, const Vector3f& out, const Vector3f& in) {
    return std::visit(evalDisneyBXDFOP{ out,in }, disneyBXDF);
}

double pdfDisneyBXDF(const DisneyBXDF& disneyBXDF, const Vector3f& out, const Vector3f& in) {
    return std::visit(pdfDisneyBXDFOP{ out,in }, disneyBXDF);
}

BSDFSampleResult sampleDisneyBXDF(const DisneyBXDF& disneyBXDF, const Vector3f& out, const Vector2f& sample) {
    BSDFSampleResult result = std::visit(sampleDisneyBXDFOP{ out,sample }, disneyBXDF);
    result.weight = evalDisneyBXDF(disneyBXDF, out, result.wi);
    result.pdf = pdfDisneyBXDF(disneyBXDF, out, result.wi);
    return result;
}

Spectrum DisneyBSDF::f(const Vector3f& out, const Vector3f& in) const {
    if (onlyDiffuse) return evalDisneyBXDF(*disneyDiffuse, out, in);
    if (onlyMetal) return evalDisneyBXDF(*disneyMetal, out, in);
    if (onlySheen) return evalDisneyBXDF(*disneySheen, out, in);
    if (onlyGlass) return evalDisneyBXDF(*disneyGlass, out, in);
    if (onlyClearCoat) return evalDisneyBXDF(*disneyClearCoat, out, in);

    double glassWeight = (1 - metallic) * specularTransmission;
    Spectrum glassResult = glassWeight > 0 ? glassWeight * evalDisneyBXDF(*disneyGlass, out, in) : Spectrum(0);
    if (CosTheta(out) < 0) {
        return glassResult;
    }

    double diffuseWeight = (1 - specularTransmission) * (1 - metallic);
    double metalWeight = (1 - specularTransmission * (1 - metallic));
    double sheenWeight = (1 - metallic) * sheen;
    double clearCoatWeight = 0.25 * clearCoat;

    Spectrum diffuseResult = diffuseWeight > 0 ? diffuseWeight * evalDisneyBXDF(*disneyDiffuse, out, in) : Spectrum(0);
    Spectrum metalResult = metalWeight > 0 ? metalWeight * evalDisneyBXDF(*disneyMetal, out, in) : Spectrum(0);
    Spectrum clearCoatResult = clearCoatWeight > 0 ? clearCoatWeight * evalDisneyBXDF(*disneyClearCoat, out, in) : Spectrum(0);
    Spectrum sheenResult = sheenWeight > 0 ? sheenWeight * evalDisneyBXDF(*disneySheen, out, in) : Spectrum(0);

    return diffuseResult + metalResult + clearCoatResult + glassResult + sheenResult;
}



BSDFSampleResult DisneyBSDF::sample(const Vector3f& out, const Vector2f& sample) const {

    if (onlyDiffuse) return sampleDisneyBXDF(*disneyDiffuse, out, sample);
    if (onlyMetal) return sampleDisneyBXDF(*disneyMetal, out, sample);
    if (onlySheen) return sampleDisneyBXDF(*disneySheen, out, sample);
    if (onlyGlass) return sampleDisneyBXDF(*disneyGlass, out, sample);
    if (onlyClearCoat) return sampleDisneyBXDF(*disneyClearCoat, out, sample);

    if (out.z() < 0) {
        return sampleDisneyBXDF(*disneyGlass, out, sample);
    }
    const double diffuseWeight = (1 - metallic) * (1 - specularTransmission);
    const double metalWeight = 1 - specularTransmission * (1 - metallic);
    const double glassWeight = (1 - metallic) * specularTransmission;
    const double clearCoatWeight = 0.25 * clearCoat;
    const double  weights[]{ diffuseWeight,metalWeight,glassWeight,clearCoatWeight };
    Distribution1D distrib(weights, 4);
    double  r = rand() % (100000 + 1) / (float)(100000 + 1); //todo modify
    int idx = distrib.SampleDiscrete(r, nullptr);
    switch (idx) {
    case 0:
        return sampleDisneyBXDF(*disneyDiffuse, out, sample);
    case 1:
        return sampleDisneyBXDF(*disneyMetal, out, sample);
    case 2:
        return sampleDisneyBXDF(*disneyGlass, out, sample);
    case 3:
        return sampleDisneyBXDF(*disneyClearCoat, out, sample);
    }
    //shod not be reached
    return BSDFSampleResult();
}


DisneyBSDF::DisneyBSDF(const Spectrum& baseColor, double specularTransmission, double metallic, double subsurface,
    double specular, double roughness, double specularTint, double anisotropic, double sheen,
    double sheenTint, double clearCoat, double clearCoatGloss, double eta
) :
    disneyDiffuse(std::make_unique<DisneyDiffuse>(DisneyDiffuse{ baseColor,roughness,subsurface })),
    disneyMetal(std::make_unique<DisneyMetal>(DisneyMetal(baseColor, roughness, anisotropic))),
    disneyClearCoat(std::make_unique<DisneyClearCoat>(DisneyClearCoat{ clearCoatGloss })),
    disneyGlass(std::make_unique<DisneyGlass>(DisneyGlass{ baseColor,roughness,anisotropic,eta })),
    disneySheen(std::make_unique<DisneySheen>(DisneySheen{ baseColor,sheenTint }))
    , specularTransmission(specularTransmission), metallic(metallic), clearCoat(clearCoat),
    sheen(sheen), ior(eta) {
}


std::shared_ptr <BSDF> DisneyMaterial::computeBSDF(const Intersection& intersect) const {
    return std::make_shared <DisneyBSDF>(
        baseColor->eval(intersect), specularTransmission->eval(intersect),
        metallic->eval(intersect), subsurface->eval(intersect), specular->eval(intersect),
        roughness->eval(intersect), specularTint->eval(intersect),
        anisotropic->eval(intersect), sheen->eval(intersect), sheenTint->eval(intersect),
        clearCoat->eval(intersect), clearCoatGloss->eval(intersect), eta);

}

DisneyMaterial::DisneyMaterial(const Json& json) : Material(json) {
    baseColor = TextureFactory::LoadTexture<Spectrum>(getChild(json, "base_color"), Spectrum(1));
    specularTransmission = TextureFactory::LoadTexture<double>(getChild(json, "specular_transmission"), 0);
    metallic = TextureFactory::LoadTexture<double>(getChild(json, "metallic"), 0);
    subsurface = TextureFactory::LoadTexture<double>(getChild(json, "subsurface"), 0);
    specular = TextureFactory::LoadTexture<double>(getChild(json, "specular"), 0);
    roughness = TextureFactory::LoadTexture<double>(getChild(json, "roughness"), 0);
    specularTint = TextureFactory::LoadTexture<double>(getChild(json, "specular_tint"), 0);
    anisotropic = TextureFactory::LoadTexture<double>(getChild(json, "anisotropic"), 0);
    sheen = TextureFactory::LoadTexture<double>(getChild(json, "sheen"), 0);
    sheenTint = TextureFactory::LoadTexture<double>(getChild(json, "sheen_tint"), 0);
    clearCoat = TextureFactory::LoadTexture<double>(getChild(json, "clear_coat"), 0);
    clearCoatGloss = TextureFactory::LoadTexture<double>(getChild(json, "clear_coat_gloss"), 0);
    eta = getOptional(json, "eta", 1.5);
    twoSideShading = false;
}


REGISTER_CLASS(DisneyMaterial, "disney")