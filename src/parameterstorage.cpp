#include "parameterstorage.h"
/*!
    - reads input file "filename" and creates parameter map
    - converts all length to internal length scale "R"
 */
ParameterStorage::ParameterStorage(std::string filename)
{
    DEBUG_FUNC_START

    std::ifstream infile(filename);
    std::string line;
    double val; // buffer for values
    std::string name; // buffer for strings

    while (std::getline(infile, line)) {
        // std::cout<<line<<std::endl;
        // remove comments
        auto pos = line.find("#");
        if (pos == 0 or line.size() == 0)
            continue;
        if (pos == std::string::npos)
            ;
        else {
            line.erase(pos, std::string::npos);
        }

        std::istringstream iss(line);
        if (!(iss >> name))
            throw std::invalid_argument("cant read line: " + line);
        if (name == "electrode") {
            double pos, voltage;
            int edge;
            if (geometry == "rect") {
                if (!(iss >> pos >> edge >> voltage))
                    throw std::invalid_argument("cant read electrode: " + line);
                electrodes.push_back({ pos, edge, voltage });
            } else if (geometry == "circle") {
                if (!(iss >> pos >> voltage))
                    throw std::invalid_argument("cant read electrode: " + line);
                electrodes.push_back({ pos, 0, voltage });
            } else {
                throw std::invalid_argument("'geometry' needs to be set first");
            }
        } else if (name == "gate") {
            if (!(iss >> gate))
                throw std::invalid_argument("cant read gate: " + line);
        } else if (name == "geometry") {
            if (!(iss >> geometry))
                throw std::invalid_argument("cant read geometry: " + line);
            if (geometry != "rect" and geometry != "circle")
                throw std::invalid_argument("invalid value for geometry: '" + line + "' valid: rect, circle");
        } else if (name == "isolateElectrode") {
            if (!(iss >> val))
                throw std::invalid_argument("cant read isolateElectrode: " + line);
            isolatedElectrodes.push_back(val);
        } else {
            if (!(iss >> val))
                throw std::invalid_argument("can't read line: " + line);
            parameters[name] = val;
            // std::cout<<name<<parameters[name]<<std::endl;
        }
    }
    if (parameters.count("randomEnergyStdDev") == 0) {
        parameters["randomEnergyStdDev"] = 0;
    }

    parameters["kT"] = parameters.at("k") * parameters.at("T");
    parameters["hoppingSiteNumber"] = parameters.at("acceptorNumber") + electrodes.size();

    // convert lens in dimensions of R
    if (geometry == "rect") {
        parameters["R"] = std::sqrt(parameters["lenX"] * parameters["lenY"] / parameters["acceptorNumber"]);
    } else if (geometry == "circle") {
        parameters["R"] = std::sqrt(M_PI * parameters.at("radius") * parameters.at("radius") / parameters["acceptorNumber"]);
    }

    parameters["minHoppingDist"] = parameters.at("minHoppingDist") / parameters.at("R");
    parameters["maxHoppingDist"] = parameters.at("maxHoppingDist") / parameters.at("R");
    parameters["maxInteractionDist"] = parameters.at("maxInteractionDist") / parameters.at("R");
    parameters["a"] = parameters.at("a") / parameters.at("R");
    parameters["minDist"] = parameters.at("minDist") / parameters.at("R");
    parameters["electrodeWidth"] = parameters.at("electrodeWidth") / parameters.at("R");
    if (geometry == "rect") {
        parameters["lenX"] = parameters.at("lenX") / parameters.at("R");
        parameters["lenY"] = parameters.at("lenY") / parameters.at("R");
    } else if (geometry == "circle") {
        parameters["radius"] = parameters.at("radius") / parameters.at("R");
    }

    double I0_joules = parameters.at("e") * parameters.at("e") / (4 * M_PI * parameters.at("eps0") * parameters.at("R") * 1e-9 * parameters.at("eps"));
    double I0_ev = I0_joules / parameters.at("e");

    // convert energies in dimensions of kT
    parameters["I0"] = I0_joules / parameters.at("kT");

    std::cout << "R = " << parameters.at("R")
              << " nm, a/R = " << parameters.at("a")
              << ", eps = " << parameters.at("eps")
              << " corresponds to I0 = " << I0_ev * 1e-3 << " meV == "
              << parameters.at("I0") << "kT" << std::endl;

    if (parameters.at("voltageScanPoints") != 1) {
        parameters["voltageScanResoultion"] = (parameters.at("voltageScanMax") - parameters.at("voltageScanMin")) / (parameters.at("voltageScanPoints") - 1);
    } else {
        parameters["voltageScanResoultion"] = 0;
    }
    for (size_t i = 0; i < parameters.at("voltageScanPoints"); i++) {
        inputVoltages.push_back(parameters.at("voltageScanMin") + i * parameters.at("voltageScanResoultion"));
    }
    std::cout << "Searching on " << parameters.at("voltageScanPoints")
              << " Voltage Points: ";
    for (size_t i = 0; i < parameters.at("voltageScanPoints"); i++) {
        std::cout << inputVoltages[i] << " ";
    }
    std::cout << std::endl;

    DEBUG_FUNC_END
}
