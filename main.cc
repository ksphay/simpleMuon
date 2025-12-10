#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4RunManager.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "FTFP_BERT.hh"

#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <cstdlib>
#include <algorithm>

#ifndef MACROS_DIR
#define MACROS_DIR "."
#endif

int main(int argc, char** argv)
{
    // ----------------------------------------
    // CLI parameters and defaults
    // ----------------------------------------
    G4String macroFile = "";       // empty => interactive + vis_muon.mac
    G4double muonEnergyGeV    = 1.0;   // default energy [GeV]

    // World extents in meters (user-specified)
    G4double xm_m = -1.0;
    G4double xp_m =  1.0;
    G4double ym_m = -1.0;
    G4double yp_m =  1.0;
    G4double zm_m = -1.0;
    G4double zp_m =  1.0;

    // Source spherical params (deg, m)
    G4double thetaDeg = 180.0;  // zenith
    G4double phiDeg   = 0.0;    // azimuth
    G4double r_m      = -1.0;   // radius (m); <0 => auto
    G4double hitPlaneHalfSize_m = 0.0; // half-size (m) for hit point shift in +/-X,Y

    // NEW: cone half-angle (deg), 0 => pencil beam
    G4double coneHalfAngleDeg = 0.0;

    // ----------------------------------------
    // Parse arguments
    // ----------------------------------------
    for (int i = 1; i < argc; ++i) {
        G4String arg = argv[i];

        if ((arg == "-E" || arg == "--energy") && (i + 1) < argc) {
            muonEnergyGeV = std::atof(argv[++i]);
        }
        else if ((arg == "--xm" || arg == "-xm") && (i + 1) < argc) {
            xm_m = std::atof(argv[++i]);
        }
        else if ((arg == "--xp" || arg == "-xp") && (i + 1) < argc) {
            xp_m = std::atof(argv[++i]);
        }
        else if ((arg == "--ym" || arg == "-ym") && (i + 1) < argc) {
            ym_m = std::atof(argv[++i]);
        }
        else if ((arg == "--yp" || arg == "-yp") && (i + 1) < argc) {
            yp_m = std::atof(argv[++i]);
        }
        else if ((arg == "--zm" || arg == "-zm") && (i + 1) < argc) {
            zm_m = std::atof(argv[++i]);
        }
        else if ((arg == "--zp" || arg == "-zp") && (i + 1) < argc) {
            zp_m = std::atof(argv[++i]);
        }
        else if ((arg == "--theta" || arg == "-th") && (i + 1) < argc) {
            thetaDeg = std::atof(argv[++i]);
        }
        else if ((arg == "--phi" || arg == "-ph") && (i + 1) < argc) {
            phiDeg = std::atof(argv[++i]);
        }
        else if (arg == "--r" && (i + 1) < argc) {
            r_m = std::atof(argv[++i]);
        }
        else if ((arg == "--cone" || arg == "-cone") && (i + 1) < argc) {
            coneHalfAngleDeg = std::atof(argv[++i]);
        }
        else if ((arg == "--hit" || arg == "-hit") && (i + 1) < argc) {
            hitPlaneHalfSize_m = std::atof(argv[++i]);
        }
        else if ((arg == "-m" || arg == "--macro") && (i + 1) < argc) {
            macroFile = argv[++i];
        }
        else if (arg[0] != '-') {
            // Positional argument: try macro if ends with .mac, else energy
            if (macroFile.empty() &&
                arg.size() >= 4 &&
                arg.substr(arg.size() - 4) == ".mac")
            {
                macroFile = arg;
            } else {
                G4double val = std::atof(arg);
                if (val > 0.0) {
                    muonEnergyGeV = val;
                } else {
                    G4cout << "[main] WARNING: positional arg \"" << arg
                           << "\" is not a positive number or .mac file, ignoring."
                           << G4endl;
                }
            }
        }
        else {
            G4cout << "[main] WARNING: unknown or incomplete argument: "
                   << arg << G4endl;
        }
    }

    // Clamp energy
    if (muonEnergyGeV <= 0.0) {
        G4cout << "[main] WARNING: non-positive energy given ("
               << muonEnergyGeV << " GeV), resetting to 1 GeV." << G4endl;
        muonEnergyGeV = 1.0;
    }

    // Ensure xm < xp, etc.
    auto ensureOrdered = [](G4double& a, G4double& b, const char* label) {
        if (b <= a) {
            G4cout << "[main] WARNING: " << label
                   << ": requested max <= min (" << a << ", " << b
                   << "), swapping." << G4endl;
            std::swap(a, b);
        }
    };
    ensureOrdered(xm_m, xp_m, "X");
    ensureOrdered(ym_m, yp_m, "Y");
    ensureOrdered(zm_m, zp_m, "Z");

    // World half sizes (meters)
    G4double halfX_m = 0.5 * (xp_m - xm_m);
    G4double halfY_m = 0.5 * (yp_m - ym_m);
    G4double halfZ_m = 0.5 * (zp_m - zm_m);
    if (halfX_m <= 0.) halfX_m = 1.0;
    if (halfY_m <= 0.) halfY_m = 1.0;
    if (halfZ_m <= 0.) halfZ_m = 1.0;

    // Convert to internal units
    G4double halfX = halfX_m * m;
    G4double halfY = halfY_m * m;
    G4double halfZ = halfZ_m * m;

    // Source radius: auto if <=0
    if (r_m <= 0.) {
        r_m = 0.8 * halfZ_m; // 80% of halfZ by default
    }

    // Clamp radius to stay inside world
    G4double minHalf_m = std::min(halfX_m, std::min(halfY_m, halfZ_m));
    G4double maxR_m    = 0.99 * minHalf_m;
    if (r_m > maxR_m) {
        G4cout << "[main] WARNING: requested radius " << r_m
               << " m exceeds world; clamping to " << maxR_m << " m."
               << G4endl;
        r_m = maxR_m;
    }
    G4double rSource = r_m * m;

    // ----------------------------------------
    // Echo configuration
    // ----------------------------------------
    G4cout << "[main] Muon energy (GeV): " << muonEnergyGeV << G4endl;
    G4cout << "[main] World spans (m): "
           << "X=[" << xm_m << ", " << xp_m << "], "
           << "Y=[" << ym_m << ", " << yp_m << "], "
           << "Z=[" << zm_m << ", " << zp_m << "]" << G4endl;
    G4cout << "[main] World half sizes (m): "
           << "X=" << halfX_m << " Y=" << halfY_m << " Z=" << halfZ_m << G4endl;
    G4cout << "[main] Source: r=" << r_m << " m, "
           << "theta=" << thetaDeg << " deg, phi=" << phiDeg << " deg" << G4endl;
    G4cout << "[main] Cone half-angle (deg): " << coneHalfAngleDeg << G4endl;
    G4cout << "[main] Hit plane half-size (m): " << hitPlaneHalfSize_m << G4endl;

    // ----------------------------------------
    // Interactive vs batch
    // ----------------------------------------
    G4bool interactive = macroFile.empty();
    G4UIExecutive* ui = nullptr;
    if (interactive) {
        ui = new G4UIExecutive(argc, argv);
    }

    // ----------------------------------------
    // Run manager
    // ----------------------------------------
#ifdef G4MULTITHREADED
    auto* runManager = new G4MTRunManager();
    runManager->SetNumberOfThreads(4);  // tune as you like
#else
    auto* runManager = new G4RunManager();
#endif

    // ----------------------------------------
    // Initialization: geometry, physics, actions
    // ----------------------------------------
    runManager->SetUserInitialization(
        new DetectorConstruction(halfX, halfY, halfZ, thetaDeg, phiDeg)
    );

    runManager->SetUserInitialization(new FTFP_BERT());

    runManager->SetUserInitialization(
        new ActionInitialization(muonEnergyGeV,
                                 thetaDeg,
                                 phiDeg,
                                 rSource,
                                 coneHalfAngleDeg,
                                 hitPlaneHalfSize_m * m)
    );

    // Visualization manager
    auto* visManager = new G4VisExecutive();
    visManager->Initialize();

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (interactive) {
        // Interactive mode: execute default vis macro
        G4String macro = G4String(MACROS_DIR) + "/vis.mac";
        G4cout << "[main] Interactive mode, executing macro: " << macro << G4endl;
        UImanager->ApplyCommand("/control/execute " + macro);

        ui->SessionStart();
        delete ui;
    } else {
        // Batch mode: execute user macro
        G4cout << "[main] Batch mode, executing macro: " << macroFile << G4endl;
        UImanager->ApplyCommand("/control/execute " + macroFile);
    }

    delete visManager;
    delete runManager;
    return 0;
}
