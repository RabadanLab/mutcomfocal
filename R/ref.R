ref<-function(ref_name) UseMethod("ref");
ref.character<-function(ref_name) {
    ref<-new.env();
    load(system.file(sprintf("refs/%s/data.RData", ref_name), package="mutcomfocal"), envir=ref);
    ref<-as.list(ref);
    class(ref)<-c("mcf_ref", class(ref));
    ref;
}
ref.default<-function (ref_name) attr(ref_name, "ref");

"ref<-"<-function(x, value) UseMethod("ref<-");
"ref<-.default"<-function (x, value) { attr(x, "ref")<-value; x }


