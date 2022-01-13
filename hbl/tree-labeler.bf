RequireVersion ("2.5.15");

LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");

LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");

labeler.analysis_description = {terms.io.info :
                            "
                            Read a tree and (optionally) a list of sequence names, and annotate branches for subsequent analyses that accept branch partitions.
                            ",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "An MSA and, optionally, a tree"
                          };


io.DisplayAnalysisBanner (labeler.analysis_description);

KeywordArgument  ("tree", "The tree to annotate (Newick format)");
SetDialogPrompt  ("Load the tree to annotate (Newick format)");

labeler.tree = trees.LoadAnnotatedTopology (TRUE);
labeler.ts = labeler.tree [^"terms.trees.newick_with_lengths"];
Topology T = labeler.ts;

KeywordArgument  ("reroot", "Reroot the tree on this node ('None' to skip rerooting)", "None");
labeler.reroot = io.PromptUserForString ("Reroot the tree a sequence from this file ('None' to skip rerooting)");

if (labeler.reroot != "None") {
    labeler.ref = alignments.ReadNucleotideDataSet ("labeler.refs", labeler.reroot);
    labeler.ref = alignments.GetSequenceNames ("labeler.refs");
    labeler.ref_dict = {};
    for (i; in; labeler.ref) {
         labeler.ref_dict[i] = 1;
    }
        
    for (labeler.reroot; in; T) {
        if (labeler.ref_dict [labeler.reroot]) {
            rerooted = RerootTree (T, labeler.reroot);
            Topology T = rerooted;
            break;
        }
    }
    
}


labeler.labels = {};

for (n; in; T) {
  if (regexp.Find (n,"^REF")) {
    labeler.labels[n] = "REFERENCE";
  } else {
    if (regexp.Find (n,"^cluster")) {
        labeler.labels[n] = "CLADE";
    }
  }
}


assert (utility.Array1D (labeler.labels) > 0, "A non-empty set of selected branches is required");

labeler.core_set = labeler.labels;

label.core = utility.Array1D (labeler.labels);

console.log ("\nSelected " + label.core + " branches to label\n");

labeler.labels * ((trees.ConjunctionLabel ("T", labeler.labels))["labels"]);

console.log ("\nLabeled " + (utility.Array1D (labeler.labels) - label.core) + " additional branches\n");

KeywordArgument ("output", "Write labeled Newick tree to");
labeler.path = io.PromptUserForFilePath ("Write labeled Newick tree to");
fprintf (labeler.path, CLEAR_FILE, tree.Annotate ("T", "relabel_and_annotate", "{}", TRUE));

KeywordArgument ("output-internal", "Write labeled Newick tree (internal branches only) to");
labeler.path_internal = io.PromptUserForFilePath ("Write labeled Newick tree (internal branches only) to");
fprintf (labeler.path_internal, CLEAR_FILE, tree.Annotate ("T", "relabel_and_annotate_internal", "{}", TRUE));


function relabel_and_annotate (node_name) {
    _label = "";
    if (labeler.labels / node_name) {
        _label = "{" + labeler.labels[node_name] + "}";
    }
    return node_name + _label;
}

function relabel_and_annotate_internal (node_name) {
    _label = "";
    if (labeler.labels / node_name && labeler.core_set / node_name == FALSE) {
        _label = "{" + labeler.labels[node_name] + "}";
    }
    return node_name + _label;
}
