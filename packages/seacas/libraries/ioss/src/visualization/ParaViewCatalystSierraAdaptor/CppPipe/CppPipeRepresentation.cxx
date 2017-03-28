
#include <vtkMultiProcessController.h>
#include <vtkSMPVRepresentationProxy.h>
#include <vtkSMPropertyHelper.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositeDataIterator.h>
#include <vtkFieldData.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkSMTransferFunctionProxy.h>
#include "CppPipeRepresentation.h"

CppPipeRepresentation::CppPipeRepresentation(vtkCPDataDescription* dataDescription,
                                             vtkSMProxy* view) : CppPipeObject()
{
  this->SetDefaults();
  this->dataDescription = dataDescription;
  this->view = view;
}

CppPipeRepresentation::CppPipeRepresentation(const Json::Value& jv,
                                             vtkCPDataDescription* dataDescription,
                                             vtkSMProxy* view) : CppPipeObject(jv)
{
  this->SetDefaults();
  this->dataDescription = dataDescription;
  this->view = view;
}

void CppPipeRepresentation::SetDefaults()
{
  this->defaults["show surfaces"] = true;
  this->defaults["show edges"] = false;
  this->defaults["show bounding box"] = false;
  this->defaults["color by blockid"] = true;
  this->defaults["show color legend"] = false;
  this->defaults["color legend position"] = Json::Value(Json::arrayValue);
  this->defaults["color legend position"].append("bottom");
  this->defaults["color legend position"].append(0.05);
  this->defaults["preset color scale"] = "Default";
  this->defaults["invert color scale"] = false;
  this->defaults["show time annotation"] = false;
  this->defaults["time annotation position"] = Json::Value(Json::arrayValue);
  this->defaults["time annotation position"].append("bottom left");
  this->defaults["time annotation position"].append(1.0);
  this->defaults["show axes"] = false;
  this->defaults["show x axis label"] = true;
  this->defaults["show y axis label"] = true;
  this->defaults["show z axis label"] = true;
  this->defaults["x axis label name"] = "X";
  this->defaults["y axis label name"] = "Y";
  this->defaults["z axis label name"] = "Z";

  std::string strColorLegendsCollection = "{\
  \"Default\":{\
    \"ColorSpace\":\"Diverging\",\
    \"NanColor\":[0.25, 0.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 0.23000000000000001, 0.29899999999999999, 0.754,\
      1.0, 0.70599999999999996, 0.016, 0.14999999999999999\
    ]\
  },\
  \"Cool to Warm\":{\
    \"ColorSpace\":\"Diverging\",\
    \"NanColor\":[0.247058823529, 0.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 0.23137254902, 0.298039215686, 0.752941176471,\
      0.5, 0.865, 0.865, 0.865,\
      1.0, 0.705882352941, 0.0156862745098, 0.149019607843\
    ]\
  },\
  \"Blue to Red Rainbow\":{\
    \"ColorSpace\":\"HSV\",\
    \"NanColor\":[0.498039215686, 0.498039215686, 0.498039215686],\
    \"RGBPoints\":[\
      0.0, 0.0, 0.0, 1.0,\
      1.0, 1.0, 0.0, 0.0\
    ]\
  },\
  \"Red to Blue Rainbow\":{\
    \"ColorSpace\":\"HSV\",\
    \"NanColor\":[0.498039215686, 0.498039215686, 0.498039215686],\
    \"RGBPoints\":[\
      0.0, 1.0, 0.0, 0.0,\
      1.0, 0.0, 0.0, 1.0\
    ]\
  },\
  \"Grayscale\":{\
    \"ColorSpace\":\"RGB\",\
    \"NanColor\":[1.0, 0.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 0.0, 0.0, 0.0,\
      1.0, 1.0, 1.0, 1.0\
    ]\
  },\
  \"X Ray\":{\
    \"ColorSpace\":\"RGB\",\
    \"NanColor\":[1.0, 0.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 1.0, 1.0, 1.0,\
      1.0, 0.0, 0.0, 0.0\
    ]\
  },\
  \"Blue to Yellow\":{\
    \"ColorSpace\":\"RGB\",\
    \"NanColor\":[1.0, 0.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 0.0392156862745, 0.0392156862745, 0.949019607843,\
      1.0, 0.949019607843, 0.949019607843, 0.0392156862745\
    ]\
  },\
  \"Black Body Radiation\":{\
    \"ColorSpace\":\"RGB\",\
    \"NanColor\":[0.0, 0.498039215686, 1.0],\
    \"RGBPoints\":[\
      0.0, 0.0, 0.0, 0.0,\
      0.4, 0.901960784314, 0.0, 0.0,\
      0.8, 0.901960784314, 0.901960784314, 0.0,\
      1.0, 1.0, 1.0, 1.0\
    ]\
  },\
  \"CIELab Blue to Red\":{\
    \"ColorSpace\":\"Lab\",\
    \"NanColor\":[1.0, 1.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 0.0, 0.6, 0.749019607843,\
      1.0, 0.76862745098, 0.466666666667, 0.341176470588\
    ]\
  },\
  \"Black, Blue and White\":{\
    \"ColorSpace\":\"RGB\",\
    \"NanColor\":[1.0, 1.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 0.0, 0.0, 0.0,\
      0.333, 0.0, 0.0, 0.501960784314,\
      0.666, 0.0, 0.501960784314, 1.0,\
      1.0, 1.0, 1.0, 1.0\
    ]\
  },\
  \"Black, Orange and White\":{\
    \"ColorSpace\":\"RGB\",\
    \"NanColor\":[1.0, 1.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 0.0, 0.0, 0.0,\
      0.333, 0.501960784314, 0.0, 0.0,\
      0.666, 1.0, 0.501960784314, 0.0,\
      1.0, 1.0, 1.0, 1.0\
    ]\
  },\
  \"Cold and Hot\":{\
    \"ColorSpace\":\"RGB\",\
    \"NanColor\":[1.0, 1.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 0.0, 1.0, 1.0,\
      0.45, 0.0, 0.0, 1.0,\
      0.5, 0.0, 0.0, 0.501960784314,\
      0.55, 1.0, 0.0, 0.0,\
      1.0, 1.0, 1.0, 0.0\
    ]\
  },\
  \"Rainbow Desaturated\":{\
    \"ColorSpace\":\"RGB\",\
    \"NanColor\":[1.0, 1.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 0.278431372549, 0.278431372549, 0.858823529412,\
      0.143, 0.0, 0.0, 0.360784313725,\
      0.285, 0.0, 1.0, 1.0,\
      0.429, 0.0, 0.501960784314, 0.0,\
      0.571, 1.0, 1.0, 0.0,\
      0.714, 1.0, 0.380392156863, 0.0,\
      0.857, 0.419607843137, 0.0, 0.0,\
      1.0, 0.878431372549, 0.301960784314, 0.301960784314\
    ]\
  },\
  \"Rainbow Blended White\":{\
    \"ColorSpace\":\"RGB\",\
    \"NanColor\":[1.0, 1.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 1.0, 1.0, 1.0,\
      0.17, 0.0, 0.0, 1.0,\
      0.34, 0.0, 1.0, 1.0,\
      0.5, 0.0, 1.0, 0.0,\
      0.67, 1.0, 1.0, 0.0,\
      0.84, 1.0, 0.0, 0.0,\
      1.0, 0.878431372549, 0.0, 1.0\
    ]\
  },\
  \"Rainbow Blended Grey\":{\
    \"ColorSpace\":\"RGB\",\
    \"NanColor\":[1.0, 1.0, 0.0],\
    \"RGBPoints\":[\
      0.0, 0.317647058824, 0.341176470588, 0.43137254902,\
      0.17, 0.0, 0.0, 1.0,\
      0.34, 0.0, 1.0, 1.0,\
      0.5, 0.0, 1.0, 0.0,\
      0.67, 1.0, 1.0, 0.0,\
      0.84, 1.0, 0.0, 0.0,\
      1.0, 0.878431372549, 0.0, 1.0\
    ]\
  }\
  }";

  Json::Reader reader;
  reader.parse(strColorLegendsCollection.c_str(),
               this->colorLegendCollection);
}

void CppPipeRepresentation::InitializeProxy(vtkSMProxy* proxy)
{
  if( proxy &&
      vtkSMPVRepresentationProxy::SafeDownCast(proxy) ) 
    {
    vtkSMPVRepresentationProxy* rp = vtkSMPVRepresentationProxy::SafeDownCast(proxy);

    bool ss = this->settings["show surfaces"].asBool() || this->defaults["show surfaces"].asBool();
    bool se = this->settings["show edges"].asBool() || this->defaults["show edges"].asBool();

    if( ss && se )
      vtkSMPropertyHelper(rp, "Representation").Set("Surface With Edges");
    else if( ss && !se)
      vtkSMPropertyHelper(rp, "Representation").Set("Surface");
    else if( !ss && se)
      vtkSMPropertyHelper(rp, "Representation").Set("Wireframe");
    else if( !ss && !se)
      vtkSMPropertyHelper(rp, "Representation").Set("Points");

    if(this->settings.isMember("color by scalar"))
      {
      vtkMultiBlockDataSet* grid = vtkMultiBlockDataSet::SafeDownCast(
        dataDescription->GetInputDescriptionByName("input")->GetGrid());

      vtkDoubleArray* localArrayBounds = vtkDoubleArray::New();
      localArrayBounds->InsertNextValue(vtkDoubleArray::GetDataTypeValueMax());
      localArrayBounds->InsertNextValue(vtkDoubleArray::GetDataTypeValueMin());
      vtkDataObject::AttributeTypes t = this->GetAttributeType(grid,
                                                               settings["color by scalar"].asString().c_str(),
                                                               localArrayBounds);


      vtkDoubleArray* globalArrayBounds = vtkDoubleArray::New();
      this->GetGlobalArrayBoundsParallel(globalArrayBounds,
                                         localArrayBounds);

      rp->SetScalarColoring(this->settings["color by scalar"].asString().c_str(), t);

      vtkSMProperty* lutProperty = rp->GetProperty("LookupTable");
      vtkSMPropertyHelper lutPropertyHelper(lutProperty);
      vtkSMProxy* lutProxy = lutPropertyHelper.GetAsProxy();


      vtkSMProperty* sofProperty = rp->GetProperty("ScalarOpacityFunction");
      vtkSMPropertyHelper sofPropertyHelper(sofProperty);
      vtkSMProxy* sofProxy = sofPropertyHelper.GetAsProxy();
      double range[2];
      range[0] = globalArrayBounds->GetValue(0);
      range[1] = globalArrayBounds->GetValue(1);

      vtkSMTransferFunctionProxy::RescaleTransferFunction(lutProxy, range);
      vtkSMTransferFunctionProxy::RescaleTransferFunction(sofProxy, range);

      this->InitializeLookupTable(lutProxy,
                                  range[0],
                                  range[1]);

      if(this->settings.isMember("show color legend"))
        rp->SetScalarBarVisibility(this->view, 
                                   this->settings["show color legend"].asBool());
      else
        rp->SetScalarBarVisibility(this->view, true);

      localArrayBounds->Delete();
      globalArrayBounds->Delete();
      }

    rp->UpdatePropertyInformation();
    rp->UpdateVTKObjects();
    }
}

void CppPipeRepresentation::InitializeLookupTable(vtkSMProxy* lut, 
                                                   double min, 
                                                   double max)
{
  std::string color_scale;
  if(this->settings.isMember("preset color scale"))
    color_scale = this->settings["preset color scale"].asString();
  else
    color_scale = this->defaults["preset color scale"].asString();

  Json::Value cs = this->colorLegendCollection[color_scale];
  vtkSMPropertyHelper(lut, "ColorSpace").Set(cs["ColorSpace"].asString().c_str());

  double nc[3];
  nc[0] = cs["NanColor"][0].asDouble();
  nc[1] = cs["NanColor"][1].asDouble();
  nc[2] = cs["NanColor"][2].asDouble();
  vtkSMPropertyHelper(lut, "NanColor").Set(nc,3);

  double* rgbp = new double[cs["RGBPoints"].size()];
  for(int i=0; i<cs["RGBPoints"].size(); i+=4)
    {
    rgbp[i] = min + cs["RGBPoints"][i].asDouble()*(max - min); 
    rgbp[i+1] = cs["RGBPoints"][i+1].asDouble();
    rgbp[i+2] = cs["RGBPoints"][i+2].asDouble();
    rgbp[i+3] = cs["RGBPoints"][i+3].asDouble();
    }

  vtkSMPropertyHelper(lut, "RGBPoints").Set(rgbp, cs["RGBPoints"].size());
  delete [] rgbp;

  lut->UpdatePropertyInformation();
}

vtkDataObject::AttributeTypes CppPipeRepresentation::GetAttributeType(vtkCompositeDataSet* cds,
                                                                       const char* array_name,
                                                                       vtkDoubleArray* localArrayBounds)
{
  vtkDataObject::AttributeTypes ret = vtkDataObject::NUMBER_OF_ATTRIBUTE_TYPES;
  vtkCompositeDataIterator* iter = cds->NewIterator();
  iter->GoToFirstItem();
  while (!iter->IsDoneWithTraversal())
    {
    vtkDataObject* block = iter->GetCurrentDataObject();
    vtkCompositeDataSet* ds = vtkCompositeDataSet::SafeDownCast(block);
    if(ds)
      {
      ret = this->GetAttributeType(ds, array_name, localArrayBounds);
      }
    else
      {
      double d[2];
      vtkDataSet* dataset = vtkDataSet::SafeDownCast(block);
      if(dataset)
        {
        if(dataset->GetCellData()->HasArray(settings["color by scalar"].asString().c_str()))
          {
          ret = vtkDataObject::CELL;
          vtkDataArray* da = dataset->GetCellData()->GetArray(settings["color by scalar"].asString().c_str());
          da->GetRange(d);
          }
        else if(dataset->GetPointData()->HasArray(settings["color by scalar"].asString().c_str()))
          {
          ret = vtkDataObject::POINT;
          vtkDataArray* da = dataset->GetPointData()->GetArray(settings["color by scalar"].asString().c_str());
          da->GetRange(d);
          }
        }
      if(d[0] < localArrayBounds->GetValue(0))
        localArrayBounds->SetValue(0, d[0]);
      if(d[1] > localArrayBounds->GetValue(1))
        localArrayBounds->SetValue(1, d[1]);
      }
    iter->GoToNextItem();
    }
  iter->Delete();
  return ret;
}

void CppPipeRepresentation::GetGlobalArrayBoundsParallel(vtkDoubleArray* globalBounds,
                                                          vtkDoubleArray* localBounds)
{
  globalBounds->SetNumberOfTuples(2);
  localBounds->SetValue(0, -localBounds->GetValue(0));

  vtkMultiProcessController::GetGlobalController()->AllReduce(localBounds, globalBounds, 0);

  globalBounds->SetValue(0, -globalBounds->GetValue(0));
}
