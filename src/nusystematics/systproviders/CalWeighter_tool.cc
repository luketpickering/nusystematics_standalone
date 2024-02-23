#include "nusystematics/systproviders/CalWeighter_tool.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "RwFramework/GSyst.h"

using namespace nusyst;
using namespace systtools;

// constructor passes up configuration object to base class for generic tool
// initialization and initialises our local copies of paramIds to unconfigured
//  flag values
CalWeighter::CalWeighter(fhicl::ParameterSet const &params)
    : IGENIESystProvider_tool(params) {

  for (int rescode = 0; rescode < 18; ++rescode) {
    for (auto const &[dial_name, gdial] :
         std::vector<std::pair<std::string, genie::rew::GSyst_t>>{
             {"MaRes", genie::rew::kXSecTwkDial_MaCCRES},
             {"MvRes", genie::rew::kXSecTwkDial_MvCCRES}}) {
      dial_infos.push_back(DialInfo{
          dial_name + "_ResCode" + std::to_string(rescode), rescode, gdial});

      // pre-allocate various parmeter vectors
      pidx_Params.push_back(kParamUnhandled<size_t>);
      CVs.push_back(0);
      Variations.emplace_back();
      ReWeightEngines.emplace_back();
    }
  }
}

SystMetaData CalWeighter::BuildSystMetaData(fhicl::ParameterSet const &ps,
                                            paramId_t firstId) {
  SystMetaData smd;
  SystParamHeader responseParam;

  SystParamHeader dial_variation_template;
  if (!ParseFhiclToolConfigurationParameter(ps, "alldial",
                                            dial_variation_template)) {
    std::cout << "esfoijewifjeowjeoiwfjewo" << std::endl;
    throw;
  }

  // loop through the four named parameters that this tool provides
  for (auto const &[prettyname, rescode, geniedial] : dial_infos) {
    SystParamHeader phdr = dial_variation_template;
    phdr.prettyName = prettyname;

    // Set up parameter definition with a standard tool configuration form
    // using helper function.
    phdr.systParamId = firstId++;

    // set any parameter-specific ParamHeader metadata here
    phdr.isSplineable = true;

    // add it to the metadata list to pass back.
    smd.push_back(phdr);
  }

  // Put any options that you want to propagate to the ParamHeaders options
  tool_options.put("verbosity_level", ps.get<int>("verbosity_level", 0));

  return smd;
}

bool CalWeighter::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {
  verbosity_level = tool_options.get<int>("verbosity_level", 0);

  // grab the pre-parsed param headers object
  SystMetaData const &md = GetSystMetaData();

  // loop through the named parameters that this tool provides, check that they
  // are configured, and grab their id in the current systmetadata and set up
  // and pre-calculations/configurations required.
  for (size_t i = 0; i < dial_infos.size(); ++i) {

    auto const &prettyname = dial_infos[i].prettyname;

    if (!HasParam(md, prettyname)) {
      if (verbosity_level > 1) {
        std::cout << "[INFO]: Don't have parameter " << prettyname
                  << " in SystMetaData. Skipping configuration." << std::endl;
      }
      continue;
    }
    // if the sytsmetadata (incoming paramheaders fhicl) has the parameter
    pidx_Params[i] = GetParamIndex(md, prettyname);

    if (verbosity_level > 1) {
      std::cout << "[INFO]: Have parameter " << prettyname
                << " in SystMetaData with ParamId: " << pidx_Params[i]
                << ". Configuring." << std::endl;
    }

    auto phdr = md[pidx_Params[i]];

    CVs[i] = phdr.centralParamValue;
    Variations[i] = phdr.paramVariations;

    for (auto v : Variations[i]) {
      // instantiate a new genie::rew::GReWeightNuXSecCCRES on the back of the
      // vector
      ReWeightEngines[i].emplace_back();
      ReWeightEngines[i].back().SetMode(
          genie::rew::GReWeightNuXSecCCRES::kModeMaMv);
      ReWeightEngines[i].back().SetSystematic(dial_infos[i].geniedial, v);
      // configure it to weight events
      ReWeightEngines[i].back().Reconfigure();
      if (verbosity_level > 2) {
        std::cout << "[LOUD]: Configured GReWeightNuXSecCCRES instance for "
                     "GENIE dial: "
                  << genie::rew::GSyst::AsString(dial_infos[i].geniedial)
                  << " at variation: " << v << std::endl;
      }
    }
  }
  // returning cleanly
  return true;
}

event_unit_response_t
CalWeighter::GetEventResponse(genie::EventRecord const &ev) {

  event_unit_response_t resp;
  SystMetaData const &md = GetSystMetaData();

  // return early if this event isn't one we provide responses for
  if (!ev.Summary()->ProcInfo().IsResonant() ||
      !ev.Summary()->ProcInfo().IsWeakCC()) {
    return resp;
  }

  // loop through and calculate weights
  for (size_t i = 0; i < dial_infos.size(); ++i) {
    // this parameter wasn't configured, nothing to do
    if (pidx_Params[i] == kParamUnhandled<size_t>) {
      continue;
    }

    if (ev.Summary()->ExclTagPtr()->Resonance() != dial_infos[i].rescode) {
      continue;
    }

    auto param_id = md[pidx_Params[i]].systParamId;
    // initialize the response array with this paramId
    resp.push_back({param_id, {}});

    // loop through variations for this parameter
    for (size_t v_it = 0; v_it < Variations[i].size(); ++v_it) {

      // put the response weight for this variation of this parameter into the
      // response object
      resp.back().responses.push_back(ReWeightEngines[i][v_it].CalcWeight(ev));
      if (verbosity_level > 3) {
        std::cout << "[DEBG]: For parameter " << md[pidx_Params[i]].prettyName
                  << " at variation[" << v_it << "] = " << Variations[i][v_it]
                  << " calculated weight: " << resp.back().responses.back()
                  << std::endl;
      }
    }
  }

  return resp;
}
