void readHiTree() {
  gROOT->Macro("hi.C");
  TFile file("output_test.root");
  file.cd("ana");
  TTree *tree = (TTree*)file.Get("hi");
  hi lyzer(tree);
  lyzer.Loop();
}
