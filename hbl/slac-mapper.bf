LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");

slac_json = io.ReadFromOrCreate ("Load SLAC fit file", {})["value"];

site_count = (slac_json["input"])["number of sites"];
tree_str = ((slac_json["input"])["trees"])[0];

Topology T = tree_str;

// store node order

branch_names = BranchName (T,-1);
branch_names_map = {};

slac_json = (slac_json["branch attributes"])[0];
site_counts = {};
bl = {};

for (i; in; branch_names) {
    branch_names_map + i;
    bl[i] = (slac_json[i])["Global MG94xREV"]; 
}

parent_index = Abs (T);

map_by_site  = {};
subs_by_site = {};


for (s = 0; s < site_count; s+=1) {
    site_map = {};
    subs = {{0,0,0,0}};
    site_subs = {};
    for (i,b; in; branch_names)    {
        p = parent_index[i];
        is_i = 2*((b[0][3] && 1)=="NODE");
        my_state = ((slac_json[b])["codon"])[0][s];
        SS = ((slac_json[b])["synonymous substitution count"])[0][s];
        NSS = ((slac_json[b])["nonsynonymous substitution count"])[0][s];
        
        subs[0+is_i] += SS;
        subs[1+is_i] += NSS;
        
        
        if (p >= 0) {
            p_state = ((slac_json[branch_names[p]])["codon"])[0][s];
            bn = b;
        } else {
            p_state = "";
            bn = "root";
        }
        if (SS + NSS) {
            site_subs [bn] = {{SS__,NSS__}};
        }
         if (my_state != p_state) {
            site_map[bn] = my_state;
        }
    }
    map_by_site[s] = site_map;
    site_counts[s] = subs;
    if (Abs (site_subs)) {
        subs_by_site[s] = site_subs;
    }
}

io.SpoolJSON ({"L" : bl, "SB2" : subs_by_site, "T" : tree_str, "map" : map_by_site, "subs" : site_counts}, io.PromptUserForFilePath("Save JSON to"));