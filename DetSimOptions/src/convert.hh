#include <fstream>
#include <string>

{
std::vector<string> ev_no = {"eV","1"};
std::vector<string> ev_m = {"eV","m"};
std::vector<string> m_m = {"m","m"};
std::vector<string> no_ns ={"1","ns"};

std::map<std::string, std::vector<string>> key_to_unit;
key_to_unit.insert(std::pair<std::string,std::vector<string>>("RINDEX",ev_no));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("ABSLENGTH",ev_m));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("FASTCOMPONENT",ev_no));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("SLOWCOMPONENT",ev_no));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("REEMISSIONPROB",ev_no));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("RAYLEIGH",ev_m));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("KINDEX",ev_no));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("EFFICIENCY",ev_no));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("REFLECTIVITY",ev_no));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("THICKNESS",m_m));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("PPOABSLENGTH",ev_m));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("PPOREEMISSIONPROB",ev_no));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("PPOCOMPONENT",ev_no));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("PPOTIMECONSTANT",no_ns));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("bisMSBABSLENGTH",ev_m));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("bisMSBREEMISSIONPROB",ev_no));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("bisMSBCOMPONENT",ev_no));
key_to_unit.insert(std::pair<std::string,std::vector<string>>("bisMSBTIMECONSTANT",no_ns));
std::cout<<"test for key ::"<<key_to_unit["RINDEX"][0]<<"   "<<key_to_unit["RINDEX"][1];

std::vector<const char *>  key={
"RINDEX" ,
"ABSLENGTH" ,
"FASTCOMPONENT" ,
"SLOWCOMPONENT" ,
"REEMISSIONPROB" ,
"RAYLEIGH" ,
"REFLECTIVITY",
 "KINDEX" ,
 "EFFICIENCY" ,
"THICKNESS" ,
"PPOABSLENGTH" ,
"PPOREEMISSIONPROB" ,
"PPOCOMPONENT" ,
"PPOTIMECONSTANT" ,
"bisMSBABSLENGTH" ,
"bisMSBREEMISSIONPROB" ,
"bisMSBCOMPONENT" ,
"bisMSBTIMECONSTANT" 
};

if(key.size() != key_to_unit.size()){
   std::cout<<"Error!"<<std::endl;
   exit(0);
}


const G4MaterialTable* theMaterialTable =  G4Material::GetMaterialTable();
G4int numOfMaterials = G4Material::GetNumberOfMaterials();
   
    for ( G4int i = 0 ; i < numOfMaterials ; i++){
          G4Material* aMaterial = (*theMaterialTable)[i];
          G4MaterialPropertiesTable* aMaterialPropertiesTable =  aMaterial->GetMaterialPropertiesTable();
          if (aMaterialPropertiesTable) {
              std::string make= "mkdir -p  ";
              std::string dir = "/tmp/huyuxiang/Material/" ;
              std::string m_name = aMaterial->GetName();
              std::string command = make + dir + m_name;
              const char * command_ = command.c_str();
              system(command_);

              for(G4int j = 0 ; j < key.size() ; j++){
                  G4MaterialPropertyVector* theVector =
                        aMaterialPropertiesTable->GetProperty(key[j]);
                 if(theVector){
                       std::cout<<" Material = " << aMaterial->GetName();
                       std::cout<<" property = " << key[j] << std::endl;
                       std::string key_s(key[j]);
                       std::string filename = dir + m_name +"/"+key_s;
                         
                       ofstream  afile;
                       afile.open(filename);
                        
                       std::vector<std::string> unit;
                       unit = key_to_unit[key_s];
                       std::cout<<"unit[0]="<<unit[0]<<"  unit[1]= " <<unit[1]<<std::endl;                      

                      // std::string unit1;
                      // std::string unit2;
                       if(unit[0]=="eV" && unit[1]=="1"){
                            for(int k = 0 ; k < theVector->GetVectorLength() ; k++){
                            afile <<std::left<<std::setw(20)<<std::setprecision(10)<<theVector->Energy(k)/eV<<std::left<<std::setw(6)<<"*eV"<<std::left<<std::setw(20)<<std::setprecision(10)<<(*theVector)[k]<< std::endl;
                            }        
                       }
                 
                      if(unit[0]=="eV" && unit[1]=="m"){
                            for(int k = 0 ; k < theVector->GetVectorLength() ; k++){
                            afile <<std::left<<std::setw(20)<<std::setprecision(10)<<theVector->Energy(k)/eV<<std::left<<std::setw(6)<<"*eV"<<std::left<<std::setw(20)<<std::setprecision(10)<<(*theVector)[k]/m<<std::left<<std::setw(6)<<"*m"<<std::endl;
                            }
                       }
                      if(unit[0]=="m" && unit[1]=="m"){
                            for(int k = 0 ; k < theVector->GetVectorLength() ; k++){
                            afile <<std::left<<std::setw(20)<<std::setprecision(10)<<theVector->Energy(k)/m<<std::left<<std::setw(6)<<"*m"<<std::left<<std::setw(20)<<std::setprecision(10)<<(*theVector)[k]/m<<std::left<<std::setw(6)<<"*m"<<std::endl;
                            }
                       }
                      if(unit[0]=="1" && unit[1]=="ns"){
                            for(int k = 0 ; k < theVector->GetVectorLength() ; k++){
                            afile <<std::left<<std::setw(20)<<std::setprecision(10)<<theVector->Energy(k)<<std::left<<std::setw(20)<<std::setprecision(10)<<(*theVector)[k]/ns<<std::left<<std::setw(6)<<"*ns"<<std::endl;
                            }
                       }

                       afile.close();
                            

                  }
                  
      
              }          

          }          

    }  


}

/*
std::string worktop = getenv("WORKTOP");
std::string dir_base =worktop + "/data/Simulation/DetSim/Material/LS/";
ofstream  afile;


dir_base += "PPOABSLENGTH";
afile.open(dir_base);
for(int i = 0 ; i < 779 ; i++){
    afile <<std::left<<std::setw(20)<<PPOABSEnergy[i]/eV<<std::left<<std::setw(6) << "*eV"<<std::left<<std::setw(20)<<PPOABSLENGTH[i]/m <<std::left<<std::setw(6) << "*m"<< std::endl;

}
afile.close();


dir_base =worktop + "/data/Simulation/DetSim/Material/LS/";
dir_base += "PPOREEMISSIONPROB";
afile.open(dir_base);
for(int i = 0 ; i < 15 ; i++){
    afile <<std::left<<std::setw(20)<<PPO_ReemMomentum[i]/eV<<std::left<<std::setw(6) << "*eV"<<std::left<<std::setw(20)<<PPO_Reemission[i] << std::endl;

}
afile.close();

dir_base =worktop + "/data/Simulation/DetSim/Material/LS/";
dir_base += "PPOCOMPONENT";
afile.open(dir_base);
for(int i = 0 ; i < 302 ; i++){
    afile <<std::left<<std::setw(20)<<PPOCOMEnergy[i]/eV<<std::left<<std::setw(6) << "*eV"<<std::left<<std::setw(20)<<PPOCOMPONENT[i] << std::endl;

}
afile.close();

dir_base =worktop + "/data/Simulation/DetSim/Material/LS/";
dir_base += "PPOCOMPONENT";
afile.open(dir_base);
for(int i = 0 ; i < 302 ; i++){
    afile <<std::left<<std::setw(20)<<PPOCOMEnergy[i]/eV<<std::left<<std::setw(6) << "*eV"<<std::left<<std::setw(20)<<PPOCOMPONENT[i] << std::endl;

}
afile.close();
*/



