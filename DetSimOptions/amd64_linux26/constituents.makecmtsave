
#-- start of constituents_header ------

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

tags      = $(tag),$(CMTEXTRATAGS)

DetSimOptions_tag = $(tag)

#cmt_local_tagfile = $(DetSimOptions_tag).make
cmt_local_tagfile = $(bin)$(DetSimOptions_tag).make

#-include $(cmt_local_tagfile)
include $(cmt_local_tagfile)

#cmt_local_setup = $(bin)setup$$$$.make
#cmt_local_setup = $(bin)$(package)setup$$$$.make
#cmt_final_setup = $(bin)DetSimOptionssetup.make
cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)$(package)setup.make

cmt_build_library_linksstamp = $(bin)cmt_build_library_links.stamp
#--------------------------------------------------------

#cmt_lock_setup = /tmp/lock$(cmt_lock_pid).make
#cmt_temp_tag = /tmp/tag$(cmt_lock_pid).make

#first :: $(cmt_local_tagfile)
#	@echo $(cmt_local_tagfile) ok
#ifndef QUICK
#first :: $(cmt_final_setup) ;
#else
#first :: ;
#endif

##	@bin=`$(cmtexe) show macro_value bin`

#$(cmt_local_tagfile) : $(cmt_lock_setup)
#	@echo "#CMT> Error: $@: No such file" >&2; exit 1
#$(cmt_local_tagfile) :
#	@echo "#CMT> Warning: $@: No such file" >&2; exit
#	@echo "#CMT> Info: $@: No need to rebuild file" >&2; exit

#$(cmt_final_setup) : $(cmt_local_tagfile) 
#	$(echo) "(constituents.make) Rebuilding $@"
#	@if test ! -d $(@D); then $(mkdir) -p $(@D); fi; \
#	  if test -f $(cmt_local_setup); then /bin/rm -f $(cmt_local_setup); fi; \
#	  trap '/bin/rm -f $(cmt_local_setup)' 0 1 2 15; \
#	  $(cmtexe) -tag=$(tags) show setup >>$(cmt_local_setup); \
#	  if test ! -f $@; then \
#	    mv $(cmt_local_setup) $@; \
#	  else \
#	    if /usr/bin/diff $(cmt_local_setup) $@ >/dev/null ; then \
#	      : ; \
#	    else \
#	      mv $(cmt_local_setup) $@; \
#	    fi; \
#	  fi

#	@/bin/echo $@ ok   

#config :: checkuses
#	@exit 0
#checkuses : ;

env.make ::
	printenv >env.make.tmp; $(cmtexe) check files env.make.tmp env.make

ifndef QUICK
all :: build_library_links ;
else
all :: $(cmt_build_library_linksstamp) ;
endif

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

dirs :: requirements
	@if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi
#	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
#	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

#requirements :
#	@if test ! -r requirements ; then echo "No requirements file" ; fi

build_library_links : dirs
	$(echo) "(constituents.make) Rebuilding library links"; \
	 $(build_library_links)
#	if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi; \
#	$(build_library_links)

$(cmt_build_library_linksstamp) : $(cmt_final_setup) $(cmt_local_tagfile) $(bin)library_links.in
	$(echo) "(constituents.make) Rebuilding library links"; \
	 $(build_library_links) -f=$(bin)library_links.in -without_cmt
	$(silent) \touch $@

ifndef PEDANTIC
.DEFAULT ::
#.DEFAULT :
	$(echo) "(constituents.make) $@: No rule for such target" >&2
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of constituents_header ------
#-- start of group ------

all_groups :: all

all :: $(all_dependencies)  $(all_pre_constituents) $(all_constituents)  $(all_post_constituents)
	$(echo) "all ok."

#	@/bin/echo " all ok."

clean :: allclean

allclean ::  $(all_constituentsclean)
	$(echo) $(all_constituentsclean)
	$(echo) "allclean ok."

#	@echo $(all_constituentsclean)
#	@/bin/echo " allclean ok."

#-- end of group ------
#-- start of group ------

all_groups :: cmt_actions

cmt_actions :: $(cmt_actions_dependencies)  $(cmt_actions_pre_constituents) $(cmt_actions_constituents)  $(cmt_actions_post_constituents)
	$(echo) "cmt_actions ok."

#	@/bin/echo " cmt_actions ok."

clean :: allclean

cmt_actionsclean ::  $(cmt_actions_constituentsclean)
	$(echo) $(cmt_actions_constituentsclean)
	$(echo) "cmt_actionsclean ok."

#	@echo $(cmt_actions_constituentsclean)
#	@/bin/echo " cmt_actionsclean ok."

#-- end of group ------
#-- start of constituent_app_lib ------

cmt_DetSimOptions_has_no_target_tag = 1
cmt_DetSimOptions_has_prototypes = 1

#--------------------------------------

ifdef cmt_DetSimOptions_has_target_tag

cmt_local_tagfile_DetSimOptions = $(bin)$(DetSimOptions_tag)_DetSimOptions.make
cmt_final_setup_DetSimOptions = $(bin)setup_DetSimOptions.make
cmt_local_DetSimOptions_makefile = $(bin)DetSimOptions.make

DetSimOptions_extratags = -tag_add=target_DetSimOptions

else

cmt_local_tagfile_DetSimOptions = $(bin)$(DetSimOptions_tag).make
cmt_final_setup_DetSimOptions = $(bin)setup.make
cmt_local_DetSimOptions_makefile = $(bin)DetSimOptions.make

endif

not_DetSimOptionscompile_dependencies = { n=0; for p in $?; do m=0; for d in $(DetSimOptionscompile_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
DetSimOptionsdirs :
	@if test ! -d $(bin)DetSimOptions; then $(mkdir) -p $(bin)DetSimOptions; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)DetSimOptions
else
DetSimOptionsdirs : ;
endif

ifdef cmt_DetSimOptions_has_target_tag

ifndef QUICK
$(cmt_local_DetSimOptions_makefile) : $(DetSimOptionscompile_dependencies) build_library_links
	$(echo) "(constituents.make) Building DetSimOptions.make"; \
	  $(cmtexe) -tag=$(tags) $(DetSimOptions_extratags) build constituent_config -out=$(cmt_local_DetSimOptions_makefile) DetSimOptions
else
$(cmt_local_DetSimOptions_makefile) : $(DetSimOptionscompile_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_DetSimOptions) ] || \
	  [ ! -f $(cmt_final_setup_DetSimOptions) ] || \
	  $(not_DetSimOptionscompile_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building DetSimOptions.make"; \
	  $(cmtexe) -tag=$(tags) $(DetSimOptions_extratags) build constituent_config -out=$(cmt_local_DetSimOptions_makefile) DetSimOptions; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_DetSimOptions_makefile) : $(DetSimOptionscompile_dependencies) build_library_links
	$(echo) "(constituents.make) Building DetSimOptions.make"; \
	  $(cmtexe) -f=$(bin)DetSimOptions.in -tag=$(tags) $(DetSimOptions_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_DetSimOptions_makefile) DetSimOptions
else
$(cmt_local_DetSimOptions_makefile) : $(DetSimOptionscompile_dependencies) $(cmt_build_library_linksstamp) $(bin)DetSimOptions.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_DetSimOptions) ] || \
	  [ ! -f $(cmt_final_setup_DetSimOptions) ] || \
	  $(not_DetSimOptionscompile_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building DetSimOptions.make"; \
	  $(cmtexe) -f=$(bin)DetSimOptions.in -tag=$(tags) $(DetSimOptions_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_DetSimOptions_makefile) DetSimOptions; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(DetSimOptions_extratags) build constituent_makefile -out=$(cmt_local_DetSimOptions_makefile) DetSimOptions

DetSimOptions :: DetSimOptionscompile DetSimOptionsinstall ;

ifdef cmt_DetSimOptions_has_prototypes

DetSimOptionsprototype : $(DetSimOptionsprototype_dependencies) $(cmt_local_DetSimOptions_makefile) dirs DetSimOptionsdirs
	$(echo) "(constituents.make) Starting $@"
	@if test -f $(cmt_local_DetSimOptions_makefile); then \
	  $(MAKE) -f $(cmt_local_DetSimOptions_makefile) $@; \
	  fi
#	@$(MAKE) -f $(cmt_local_DetSimOptions_makefile) $@
	$(echo) "(constituents.make) $@ done"

DetSimOptionscompile : DetSimOptionsprototype

endif

DetSimOptionscompile : $(DetSimOptionscompile_dependencies) $(cmt_local_DetSimOptions_makefile) dirs DetSimOptionsdirs
	$(echo) "(constituents.make) Starting $@"
	@if test -f $(cmt_local_DetSimOptions_makefile); then \
	  $(MAKE) -f $(cmt_local_DetSimOptions_makefile) $@; \
	  fi
#	@$(MAKE) -f $(cmt_local_DetSimOptions_makefile) $@
	$(echo) "(constituents.make) $@ done"

clean :: DetSimOptionsclean ;

DetSimOptionsclean :: $(DetSimOptionsclean_dependencies) ##$(cmt_local_DetSimOptions_makefile)
	$(echo) "(constituents.make) Starting $@"
	@-if test -f $(cmt_local_DetSimOptions_makefile); then \
	  $(MAKE) -f $(cmt_local_DetSimOptions_makefile) $@; \
	fi
	$(echo) "(constituents.make) $@ done"
#	@-$(MAKE) -f $(cmt_local_DetSimOptions_makefile) DetSimOptionsclean

##	  /bin/rm -f $(cmt_local_DetSimOptions_makefile) $(bin)DetSimOptions_dependencies.make

install :: DetSimOptionsinstall ;

DetSimOptionsinstall :: DetSimOptionscompile $(DetSimOptions_dependencies) $(cmt_local_DetSimOptions_makefile)
	$(echo) "(constituents.make) Starting $@"
	@if test -f $(cmt_local_DetSimOptions_makefile); then \
	  $(MAKE) -f $(cmt_local_DetSimOptions_makefile) $@; \
	  fi
#	@$(MAKE) -f $(cmt_local_DetSimOptions_makefile) $@
	$(echo) "(constituents.make) $@ done"

uninstall : DetSimOptionsuninstall

$(foreach d,$(DetSimOptions_dependencies),$(eval $(d)uninstall_dependencies += DetSimOptionsuninstall))

DetSimOptionsuninstall : $(DetSimOptionsuninstall_dependencies) ##$(cmt_local_DetSimOptions_makefile)
	$(echo) "(constituents.make) Starting $@"
	@-if test -f $(cmt_local_DetSimOptions_makefile); then \
	  $(MAKE) -f $(cmt_local_DetSimOptions_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_DetSimOptions_makefile) uninstall
	$(echo) "(constituents.make) $@ done"

remove_library_links :: DetSimOptionsuninstall ;

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ DetSimOptions"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ DetSimOptions done"
endif

#-- end of constituent_app_lib ------
#-- start of constituent ------

cmt_DetSimOptions_python_has_no_target_tag = 1

#--------------------------------------

ifdef cmt_DetSimOptions_python_has_target_tag

cmt_local_tagfile_DetSimOptions_python = $(bin)$(DetSimOptions_tag)_DetSimOptions_python.make
cmt_final_setup_DetSimOptions_python = $(bin)setup_DetSimOptions_python.make
cmt_local_DetSimOptions_python_makefile = $(bin)DetSimOptions_python.make

DetSimOptions_python_extratags = -tag_add=target_DetSimOptions_python

else

cmt_local_tagfile_DetSimOptions_python = $(bin)$(DetSimOptions_tag).make
cmt_final_setup_DetSimOptions_python = $(bin)setup.make
cmt_local_DetSimOptions_python_makefile = $(bin)DetSimOptions_python.make

endif

not_DetSimOptions_python_dependencies = { n=0; for p in $?; do m=0; for d in $(DetSimOptions_python_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
DetSimOptions_pythondirs :
	@if test ! -d $(bin)DetSimOptions_python; then $(mkdir) -p $(bin)DetSimOptions_python; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)DetSimOptions_python
else
DetSimOptions_pythondirs : ;
endif

ifdef cmt_DetSimOptions_python_has_target_tag

ifndef QUICK
$(cmt_local_DetSimOptions_python_makefile) : $(DetSimOptions_python_dependencies) build_library_links
	$(echo) "(constituents.make) Building DetSimOptions_python.make"; \
	  $(cmtexe) -tag=$(tags) $(DetSimOptions_python_extratags) build constituent_config -out=$(cmt_local_DetSimOptions_python_makefile) DetSimOptions_python
else
$(cmt_local_DetSimOptions_python_makefile) : $(DetSimOptions_python_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_DetSimOptions_python) ] || \
	  [ ! -f $(cmt_final_setup_DetSimOptions_python) ] || \
	  $(not_DetSimOptions_python_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building DetSimOptions_python.make"; \
	  $(cmtexe) -tag=$(tags) $(DetSimOptions_python_extratags) build constituent_config -out=$(cmt_local_DetSimOptions_python_makefile) DetSimOptions_python; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_DetSimOptions_python_makefile) : $(DetSimOptions_python_dependencies) build_library_links
	$(echo) "(constituents.make) Building DetSimOptions_python.make"; \
	  $(cmtexe) -f=$(bin)DetSimOptions_python.in -tag=$(tags) $(DetSimOptions_python_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_DetSimOptions_python_makefile) DetSimOptions_python
else
$(cmt_local_DetSimOptions_python_makefile) : $(DetSimOptions_python_dependencies) $(cmt_build_library_linksstamp) $(bin)DetSimOptions_python.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_DetSimOptions_python) ] || \
	  [ ! -f $(cmt_final_setup_DetSimOptions_python) ] || \
	  $(not_DetSimOptions_python_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building DetSimOptions_python.make"; \
	  $(cmtexe) -f=$(bin)DetSimOptions_python.in -tag=$(tags) $(DetSimOptions_python_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_DetSimOptions_python_makefile) DetSimOptions_python; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(DetSimOptions_python_extratags) build constituent_makefile -out=$(cmt_local_DetSimOptions_python_makefile) DetSimOptions_python

DetSimOptions_python :: $(DetSimOptions_python_dependencies) $(cmt_local_DetSimOptions_python_makefile) dirs DetSimOptions_pythondirs
	$(echo) "(constituents.make) Starting DetSimOptions_python"
	@if test -f $(cmt_local_DetSimOptions_python_makefile); then \
	  $(MAKE) -f $(cmt_local_DetSimOptions_python_makefile) DetSimOptions_python; \
	  fi
#	@$(MAKE) -f $(cmt_local_DetSimOptions_python_makefile) DetSimOptions_python
	$(echo) "(constituents.make) DetSimOptions_python done"

clean :: DetSimOptions_pythonclean ;

DetSimOptions_pythonclean :: $(DetSimOptions_pythonclean_dependencies) ##$(cmt_local_DetSimOptions_python_makefile)
	$(echo) "(constituents.make) Starting DetSimOptions_pythonclean"
	@-if test -f $(cmt_local_DetSimOptions_python_makefile); then \
	  $(MAKE) -f $(cmt_local_DetSimOptions_python_makefile) DetSimOptions_pythonclean; \
	fi
	$(echo) "(constituents.make) DetSimOptions_pythonclean done"
#	@-$(MAKE) -f $(cmt_local_DetSimOptions_python_makefile) DetSimOptions_pythonclean

##	  /bin/rm -f $(cmt_local_DetSimOptions_python_makefile) $(bin)DetSimOptions_python_dependencies.make

install :: DetSimOptions_pythoninstall ;

DetSimOptions_pythoninstall :: $(DetSimOptions_python_dependencies) $(cmt_local_DetSimOptions_python_makefile)
	$(echo) "(constituents.make) Starting $@"
	@if test -f $(cmt_local_DetSimOptions_python_makefile); then \
	  $(MAKE) -f $(cmt_local_DetSimOptions_python_makefile) install; \
	  fi
#	@-$(MAKE) -f $(cmt_local_DetSimOptions_python_makefile) install
	$(echo) "(constituents.make) $@ done"

uninstall : DetSimOptions_pythonuninstall

$(foreach d,$(DetSimOptions_python_dependencies),$(eval $(d)uninstall_dependencies += DetSimOptions_pythonuninstall))

DetSimOptions_pythonuninstall : $(DetSimOptions_pythonuninstall_dependencies) ##$(cmt_local_DetSimOptions_python_makefile)
	$(echo) "(constituents.make) Starting $@"
	@-if test -f $(cmt_local_DetSimOptions_python_makefile); then \
	  $(MAKE) -f $(cmt_local_DetSimOptions_python_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_DetSimOptions_python_makefile) uninstall
	$(echo) "(constituents.make) $@ done"

remove_library_links :: DetSimOptions_pythonuninstall ;

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ DetSimOptions_python"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ DetSimOptions_python done"
endif

#-- end of constituent ------
#-- start of constituent ------

cmt_make_has_no_target_tag = 1

#--------------------------------------

ifdef cmt_make_has_target_tag

cmt_local_tagfile_make = $(bin)$(DetSimOptions_tag)_make.make
cmt_final_setup_make = $(bin)setup_make.make
cmt_local_make_makefile = $(bin)make.make

make_extratags = -tag_add=target_make

else

cmt_local_tagfile_make = $(bin)$(DetSimOptions_tag).make
cmt_final_setup_make = $(bin)setup.make
cmt_local_make_makefile = $(bin)make.make

endif

not_make_dependencies = { n=0; for p in $?; do m=0; for d in $(make_dependencies); do if [ $$p = $$d ]; then m=1; break; fi; done; if [ $$m -eq 0 ]; then n=1; break; fi; done; [ $$n -eq 1 ]; }

ifdef STRUCTURED_OUTPUT
makedirs :
	@if test ! -d $(bin)make; then $(mkdir) -p $(bin)make; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)make
else
makedirs : ;
endif

ifdef cmt_make_has_target_tag

ifndef QUICK
$(cmt_local_make_makefile) : $(make_dependencies) build_library_links
	$(echo) "(constituents.make) Building make.make"; \
	  $(cmtexe) -tag=$(tags) $(make_extratags) build constituent_config -out=$(cmt_local_make_makefile) make
else
$(cmt_local_make_makefile) : $(make_dependencies) $(cmt_build_library_linksstamp) $(use_requirements)
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_make) ] || \
	  [ ! -f $(cmt_final_setup_make) ] || \
	  $(not_make_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building make.make"; \
	  $(cmtexe) -tag=$(tags) $(make_extratags) build constituent_config -out=$(cmt_local_make_makefile) make; \
	  fi
endif

else

ifndef QUICK
$(cmt_local_make_makefile) : $(make_dependencies) build_library_links
	$(echo) "(constituents.make) Building make.make"; \
	  $(cmtexe) -f=$(bin)make.in -tag=$(tags) $(make_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_make_makefile) make
else
$(cmt_local_make_makefile) : $(make_dependencies) $(cmt_build_library_linksstamp) $(bin)make.in
	@if [ ! -f $@ ] || [ ! -f $(cmt_local_tagfile_make) ] || \
	  [ ! -f $(cmt_final_setup_make) ] || \
	  $(not_make_dependencies) ; then \
	  test -z "$(cmtmsg)" || \
	  echo "$(CMTMSGPREFIX)" "(constituents.make) Building make.make"; \
	  $(cmtexe) -f=$(bin)make.in -tag=$(tags) $(make_extratags) build constituent_makefile -without_cmt -out=$(cmt_local_make_makefile) make; \
	  fi
endif

endif

#	  $(cmtexe) -tag=$(tags) $(make_extratags) build constituent_makefile -out=$(cmt_local_make_makefile) make

make :: $(make_dependencies) $(cmt_local_make_makefile) dirs makedirs
	$(echo) "(constituents.make) Starting make"
	@if test -f $(cmt_local_make_makefile); then \
	  $(MAKE) -f $(cmt_local_make_makefile) make; \
	  fi
#	@$(MAKE) -f $(cmt_local_make_makefile) make
	$(echo) "(constituents.make) make done"

clean :: makeclean ;

makeclean :: $(makeclean_dependencies) ##$(cmt_local_make_makefile)
	$(echo) "(constituents.make) Starting makeclean"
	@-if test -f $(cmt_local_make_makefile); then \
	  $(MAKE) -f $(cmt_local_make_makefile) makeclean; \
	fi
	$(echo) "(constituents.make) makeclean done"
#	@-$(MAKE) -f $(cmt_local_make_makefile) makeclean

##	  /bin/rm -f $(cmt_local_make_makefile) $(bin)make_dependencies.make

install :: makeinstall ;

makeinstall :: $(make_dependencies) $(cmt_local_make_makefile)
	$(echo) "(constituents.make) Starting $@"
	@if test -f $(cmt_local_make_makefile); then \
	  $(MAKE) -f $(cmt_local_make_makefile) install; \
	  fi
#	@-$(MAKE) -f $(cmt_local_make_makefile) install
	$(echo) "(constituents.make) $@ done"

uninstall : makeuninstall

$(foreach d,$(make_dependencies),$(eval $(d)uninstall_dependencies += makeuninstall))

makeuninstall : $(makeuninstall_dependencies) ##$(cmt_local_make_makefile)
	$(echo) "(constituents.make) Starting $@"
	@-if test -f $(cmt_local_make_makefile); then \
	  $(MAKE) -f $(cmt_local_make_makefile) uninstall; \
	  fi
#	@$(MAKE) -f $(cmt_local_make_makefile) uninstall
	$(echo) "(constituents.make) $@ done"

remove_library_links :: makeuninstall ;

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ make"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ make done"
endif

#-- end of constituent ------
#-- start of constituents_trailer ------

uninstall : remove_library_links ;
clean ::
	$(cleanup_echo) $(cmt_build_library_linksstamp)
	-$(cleanup_silent) \rm -f $(cmt_build_library_linksstamp)
#clean :: remove_library_links

remove_library_links ::
ifndef QUICK
	$(echo) "(constituents.make) Removing library links"; \
	  $(remove_library_links)
else
	$(echo) "(constituents.make) Removing library links"; \
	  $(remove_library_links) -f=$(bin)library_links.in -without_cmt
endif

#-- end of constituents_trailer ------
