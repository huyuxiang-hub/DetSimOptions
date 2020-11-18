#-- start of make_header -----------------

#====================================
#  Library DetSimOptions
#
#   Generated Mon Nov 16 16:25:36 2020  by huyuxiang
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_DetSimOptions_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_DetSimOptions_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_DetSimOptions

DetSimOptions_tag = $(tag)

#cmt_local_tagfile_DetSimOptions = $(DetSimOptions_tag)_DetSimOptions.make
cmt_local_tagfile_DetSimOptions = $(bin)$(DetSimOptions_tag)_DetSimOptions.make

else

tags      = $(tag),$(CMTEXTRATAGS)

DetSimOptions_tag = $(tag)

#cmt_local_tagfile_DetSimOptions = $(DetSimOptions_tag).make
cmt_local_tagfile_DetSimOptions = $(bin)$(DetSimOptions_tag).make

endif

include $(cmt_local_tagfile_DetSimOptions)
#-include $(cmt_local_tagfile_DetSimOptions)

ifdef cmt_DetSimOptions_has_target_tag

cmt_final_setup_DetSimOptions = $(bin)setup_DetSimOptions.make
cmt_dependencies_in_DetSimOptions = $(bin)dependencies_DetSimOptions.in
#cmt_final_setup_DetSimOptions = $(bin)DetSimOptions_DetSimOptionssetup.make
cmt_local_DetSimOptions_makefile = $(bin)DetSimOptions.make

else

cmt_final_setup_DetSimOptions = $(bin)setup.make
cmt_dependencies_in_DetSimOptions = $(bin)dependencies.in
#cmt_final_setup_DetSimOptions = $(bin)DetSimOptionssetup.make
cmt_local_DetSimOptions_makefile = $(bin)DetSimOptions.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)DetSimOptionssetup.make

#DetSimOptions :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'DetSimOptions'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = DetSimOptions/
#DetSimOptions::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

DetSimOptionslibname   = $(bin)$(library_prefix)DetSimOptions$(library_suffix)
DetSimOptionslib       = $(DetSimOptionslibname).a
DetSimOptionsstamp     = $(bin)DetSimOptions.stamp
DetSimOptionsshstamp   = $(bin)DetSimOptions.shstamp

DetSimOptions :: dirs  DetSimOptionsLIB
	$(echo) "DetSimOptions ok"

cmt_DetSimOptions_has_prototypes = 1

#--------------------------------------

ifdef cmt_DetSimOptions_has_prototypes

DetSimOptionsprototype :  ;

endif

DetSimOptionscompile : $(bin)LSExpPrimaryGeneratorAction.o $(bin)DetSim0Svc.o $(bin)WaterPoolConstruction.o $(bin)LSExpTrackingAction.o $(bin)LSExpEventAction.o $(bin)LSExpRunAction.o $(bin)ExpHallConstruction.o $(bin)RockConstruction.o $(bin)LSExpSteppingAction.o $(bin)LSExpDetectorConstruction.o $(bin)LSExpParticleGun.o $(bin)LSExpPhysicsList.o $(bin)DetPrototypeSvc.o $(bin)GlobalGeomInfo.o ;

#-- end of libary_header ----------------
#-- start of libary ----------------------

DetSimOptionsLIB :: $(DetSimOptionslib) $(DetSimOptionsshstamp)
	$(echo) "DetSimOptions : library ok"

$(DetSimOptionslib) :: $(bin)LSExpPrimaryGeneratorAction.o $(bin)DetSim0Svc.o $(bin)WaterPoolConstruction.o $(bin)LSExpTrackingAction.o $(bin)LSExpEventAction.o $(bin)LSExpRunAction.o $(bin)ExpHallConstruction.o $(bin)RockConstruction.o $(bin)LSExpSteppingAction.o $(bin)LSExpDetectorConstruction.o $(bin)LSExpParticleGun.o $(bin)LSExpPhysicsList.o $(bin)DetPrototypeSvc.o $(bin)GlobalGeomInfo.o
	$(lib_echo) "static library $@"
	$(lib_silent) [ ! -f $@ ] || \rm -f $@
	$(lib_silent) $(ar) $(DetSimOptionslib) $(bin)LSExpPrimaryGeneratorAction.o $(bin)DetSim0Svc.o $(bin)WaterPoolConstruction.o $(bin)LSExpTrackingAction.o $(bin)LSExpEventAction.o $(bin)LSExpRunAction.o $(bin)ExpHallConstruction.o $(bin)RockConstruction.o $(bin)LSExpSteppingAction.o $(bin)LSExpDetectorConstruction.o $(bin)LSExpParticleGun.o $(bin)LSExpPhysicsList.o $(bin)DetPrototypeSvc.o $(bin)GlobalGeomInfo.o
	$(lib_silent) $(ranlib) $(DetSimOptionslib)
	$(lib_silent) cat /dev/null >$(DetSimOptionsstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

#
# We add one level of dependency upon the true shared library 
# (rather than simply upon the stamp file)
# this is for cases where the shared library has not been built
# while the stamp was created (error??) 
#

$(DetSimOptionslibname).$(shlibsuffix) :: $(DetSimOptionslib) requirements $(use_requirements) $(DetSimOptionsstamps)
	$(lib_echo) "shared library $@"
	$(lib_silent) if test "$(makecmd)"; then QUIET=; else QUIET=1; fi; QUIET=$${QUIET} bin="$(bin)" ld="$(shlibbuilder)" ldflags="$(shlibflags)" suffix=$(shlibsuffix) libprefix=$(library_prefix) libsuffix=$(library_suffix) $(make_shlib) "$(tags)" DetSimOptions $(DetSimOptions_shlibflags)
	$(lib_silent) cat /dev/null >$(DetSimOptionsshstamp)

$(DetSimOptionsshstamp) :: $(DetSimOptionslibname).$(shlibsuffix)
	$(lib_silent) if test -f $(DetSimOptionslibname).$(shlibsuffix) ; then cat /dev/null >$(DetSimOptionsshstamp) ; fi

DetSimOptionsclean ::
	$(cleanup_echo) objects DetSimOptions
	$(cleanup_silent) /bin/rm -f $(bin)LSExpPrimaryGeneratorAction.o $(bin)DetSim0Svc.o $(bin)WaterPoolConstruction.o $(bin)LSExpTrackingAction.o $(bin)LSExpEventAction.o $(bin)LSExpRunAction.o $(bin)ExpHallConstruction.o $(bin)RockConstruction.o $(bin)LSExpSteppingAction.o $(bin)LSExpDetectorConstruction.o $(bin)LSExpParticleGun.o $(bin)LSExpPhysicsList.o $(bin)DetPrototypeSvc.o $(bin)GlobalGeomInfo.o
	$(cleanup_silent) /bin/rm -f $(patsubst %.o,%.d,$(bin)LSExpPrimaryGeneratorAction.o $(bin)DetSim0Svc.o $(bin)WaterPoolConstruction.o $(bin)LSExpTrackingAction.o $(bin)LSExpEventAction.o $(bin)LSExpRunAction.o $(bin)ExpHallConstruction.o $(bin)RockConstruction.o $(bin)LSExpSteppingAction.o $(bin)LSExpDetectorConstruction.o $(bin)LSExpParticleGun.o $(bin)LSExpPhysicsList.o $(bin)DetPrototypeSvc.o $(bin)GlobalGeomInfo.o) $(patsubst %.o,%.dep,$(bin)LSExpPrimaryGeneratorAction.o $(bin)DetSim0Svc.o $(bin)WaterPoolConstruction.o $(bin)LSExpTrackingAction.o $(bin)LSExpEventAction.o $(bin)LSExpRunAction.o $(bin)ExpHallConstruction.o $(bin)RockConstruction.o $(bin)LSExpSteppingAction.o $(bin)LSExpDetectorConstruction.o $(bin)LSExpParticleGun.o $(bin)LSExpPhysicsList.o $(bin)DetPrototypeSvc.o $(bin)GlobalGeomInfo.o) $(patsubst %.o,%.d.stamp,$(bin)LSExpPrimaryGeneratorAction.o $(bin)DetSim0Svc.o $(bin)WaterPoolConstruction.o $(bin)LSExpTrackingAction.o $(bin)LSExpEventAction.o $(bin)LSExpRunAction.o $(bin)ExpHallConstruction.o $(bin)RockConstruction.o $(bin)LSExpSteppingAction.o $(bin)LSExpDetectorConstruction.o $(bin)LSExpParticleGun.o $(bin)LSExpPhysicsList.o $(bin)DetPrototypeSvc.o $(bin)GlobalGeomInfo.o)
	$(cleanup_silent) cd $(bin); /bin/rm -rf DetSimOptions_deps DetSimOptions_dependencies.make

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

install_dir = ${CMTINSTALLAREA}/$(tag)/lib
DetSimOptionsinstallname = $(library_prefix)DetSimOptions$(library_suffix).$(shlibsuffix)

DetSimOptions :: DetSimOptionsinstall ;

install :: DetSimOptionsinstall ;

DetSimOptionsinstall :: $(install_dir)/$(DetSimOptionsinstallname)
ifdef CMTINSTALLAREA
	$(echo) "installation done"
endif

$(install_dir)/$(DetSimOptionsinstallname) :: $(bin)$(DetSimOptionsinstallname)
ifdef CMTINSTALLAREA
	$(install_silent) $(cmt_install_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(DetSimOptionsinstallname)" \
	    -out "$(install_dir)" \
	    -cmd "$(cmt_installarea_command)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

##DetSimOptionsclean :: DetSimOptionsuninstall

uninstall :: DetSimOptionsuninstall ;

DetSimOptionsuninstall ::
ifdef CMTINSTALLAREA
	$(cleanup_silent) $(cmt_uninstall_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(DetSimOptionsinstallname)" \
	    -out "$(install_dir)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

#-- end of libary -----------------------
#-- start of dependencies ------------------
ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
ifneq ($(MAKECMDGOALS),DetSimOptionsprototype)

$(bin)DetSimOptions_dependencies.make : $(use_requirements) $(cmt_final_setup_DetSimOptions)
	$(echo) "(DetSimOptions.make) Rebuilding $@"; \
	  $(build_dependencies) -out=$@ -start_all $(src)LSExpPrimaryGeneratorAction.cc $(src)DetSim0Svc.cc $(src)WaterPoolConstruction.cc $(src)LSExpTrackingAction.cc $(src)LSExpEventAction.cc $(src)LSExpRunAction.cc $(src)ExpHallConstruction.cc $(src)RockConstruction.cc $(src)LSExpSteppingAction.cc $(src)LSExpDetectorConstruction.cc $(src)LSExpParticleGun.cc $(src)LSExpPhysicsList.cc $(src)DetPrototypeSvc.cc $(src)GlobalGeomInfo.cc -end_all $(includes) $(app_DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) -name=DetSimOptions $? -f=$(cmt_dependencies_in_DetSimOptions) -without_cmt

-include $(bin)DetSimOptions_dependencies.make

endif
endif
endif

DetSimOptionsclean ::
	$(cleanup_silent) \rm -rf $(bin)DetSimOptions_deps $(bin)DetSimOptions_dependencies.make
#-- end of dependencies -------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)LSExpPrimaryGeneratorAction.d

$(bin)$(binobj)LSExpPrimaryGeneratorAction.d :

$(bin)$(binobj)LSExpPrimaryGeneratorAction.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)LSExpPrimaryGeneratorAction.o : $(src)LSExpPrimaryGeneratorAction.cc
	$(cpp_echo) $(src)LSExpPrimaryGeneratorAction.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpPrimaryGeneratorAction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpPrimaryGeneratorAction_cppflags) $(LSExpPrimaryGeneratorAction_cc_cppflags)  $(src)LSExpPrimaryGeneratorAction.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(LSExpPrimaryGeneratorAction_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)LSExpPrimaryGeneratorAction.cc

$(bin)$(binobj)LSExpPrimaryGeneratorAction.o : $(LSExpPrimaryGeneratorAction_cc_dependencies)
	$(cpp_echo) $(src)LSExpPrimaryGeneratorAction.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpPrimaryGeneratorAction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpPrimaryGeneratorAction_cppflags) $(LSExpPrimaryGeneratorAction_cc_cppflags)  $(src)LSExpPrimaryGeneratorAction.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)DetSim0Svc.d

$(bin)$(binobj)DetSim0Svc.d :

$(bin)$(binobj)DetSim0Svc.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)DetSim0Svc.o : $(src)DetSim0Svc.cc
	$(cpp_echo) $(src)DetSim0Svc.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(DetSim0Svc_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(DetSim0Svc_cppflags) $(DetSim0Svc_cc_cppflags)  $(src)DetSim0Svc.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(DetSim0Svc_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)DetSim0Svc.cc

$(bin)$(binobj)DetSim0Svc.o : $(DetSim0Svc_cc_dependencies)
	$(cpp_echo) $(src)DetSim0Svc.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(DetSim0Svc_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(DetSim0Svc_cppflags) $(DetSim0Svc_cc_cppflags)  $(src)DetSim0Svc.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)WaterPoolConstruction.d

$(bin)$(binobj)WaterPoolConstruction.d :

$(bin)$(binobj)WaterPoolConstruction.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)WaterPoolConstruction.o : $(src)WaterPoolConstruction.cc
	$(cpp_echo) $(src)WaterPoolConstruction.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(WaterPoolConstruction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(WaterPoolConstruction_cppflags) $(WaterPoolConstruction_cc_cppflags)  $(src)WaterPoolConstruction.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(WaterPoolConstruction_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)WaterPoolConstruction.cc

$(bin)$(binobj)WaterPoolConstruction.o : $(WaterPoolConstruction_cc_dependencies)
	$(cpp_echo) $(src)WaterPoolConstruction.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(WaterPoolConstruction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(WaterPoolConstruction_cppflags) $(WaterPoolConstruction_cc_cppflags)  $(src)WaterPoolConstruction.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)LSExpTrackingAction.d

$(bin)$(binobj)LSExpTrackingAction.d :

$(bin)$(binobj)LSExpTrackingAction.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)LSExpTrackingAction.o : $(src)LSExpTrackingAction.cc
	$(cpp_echo) $(src)LSExpTrackingAction.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpTrackingAction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpTrackingAction_cppflags) $(LSExpTrackingAction_cc_cppflags)  $(src)LSExpTrackingAction.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(LSExpTrackingAction_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)LSExpTrackingAction.cc

$(bin)$(binobj)LSExpTrackingAction.o : $(LSExpTrackingAction_cc_dependencies)
	$(cpp_echo) $(src)LSExpTrackingAction.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpTrackingAction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpTrackingAction_cppflags) $(LSExpTrackingAction_cc_cppflags)  $(src)LSExpTrackingAction.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)LSExpEventAction.d

$(bin)$(binobj)LSExpEventAction.d :

$(bin)$(binobj)LSExpEventAction.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)LSExpEventAction.o : $(src)LSExpEventAction.cc
	$(cpp_echo) $(src)LSExpEventAction.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpEventAction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpEventAction_cppflags) $(LSExpEventAction_cc_cppflags)  $(src)LSExpEventAction.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(LSExpEventAction_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)LSExpEventAction.cc

$(bin)$(binobj)LSExpEventAction.o : $(LSExpEventAction_cc_dependencies)
	$(cpp_echo) $(src)LSExpEventAction.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpEventAction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpEventAction_cppflags) $(LSExpEventAction_cc_cppflags)  $(src)LSExpEventAction.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)LSExpRunAction.d

$(bin)$(binobj)LSExpRunAction.d :

$(bin)$(binobj)LSExpRunAction.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)LSExpRunAction.o : $(src)LSExpRunAction.cc
	$(cpp_echo) $(src)LSExpRunAction.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpRunAction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpRunAction_cppflags) $(LSExpRunAction_cc_cppflags)  $(src)LSExpRunAction.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(LSExpRunAction_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)LSExpRunAction.cc

$(bin)$(binobj)LSExpRunAction.o : $(LSExpRunAction_cc_dependencies)
	$(cpp_echo) $(src)LSExpRunAction.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpRunAction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpRunAction_cppflags) $(LSExpRunAction_cc_cppflags)  $(src)LSExpRunAction.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)ExpHallConstruction.d

$(bin)$(binobj)ExpHallConstruction.d :

$(bin)$(binobj)ExpHallConstruction.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)ExpHallConstruction.o : $(src)ExpHallConstruction.cc
	$(cpp_echo) $(src)ExpHallConstruction.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(ExpHallConstruction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(ExpHallConstruction_cppflags) $(ExpHallConstruction_cc_cppflags)  $(src)ExpHallConstruction.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(ExpHallConstruction_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)ExpHallConstruction.cc

$(bin)$(binobj)ExpHallConstruction.o : $(ExpHallConstruction_cc_dependencies)
	$(cpp_echo) $(src)ExpHallConstruction.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(ExpHallConstruction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(ExpHallConstruction_cppflags) $(ExpHallConstruction_cc_cppflags)  $(src)ExpHallConstruction.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)RockConstruction.d

$(bin)$(binobj)RockConstruction.d :

$(bin)$(binobj)RockConstruction.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)RockConstruction.o : $(src)RockConstruction.cc
	$(cpp_echo) $(src)RockConstruction.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(RockConstruction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(RockConstruction_cppflags) $(RockConstruction_cc_cppflags)  $(src)RockConstruction.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(RockConstruction_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)RockConstruction.cc

$(bin)$(binobj)RockConstruction.o : $(RockConstruction_cc_dependencies)
	$(cpp_echo) $(src)RockConstruction.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(RockConstruction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(RockConstruction_cppflags) $(RockConstruction_cc_cppflags)  $(src)RockConstruction.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)LSExpSteppingAction.d

$(bin)$(binobj)LSExpSteppingAction.d :

$(bin)$(binobj)LSExpSteppingAction.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)LSExpSteppingAction.o : $(src)LSExpSteppingAction.cc
	$(cpp_echo) $(src)LSExpSteppingAction.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpSteppingAction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpSteppingAction_cppflags) $(LSExpSteppingAction_cc_cppflags)  $(src)LSExpSteppingAction.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(LSExpSteppingAction_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)LSExpSteppingAction.cc

$(bin)$(binobj)LSExpSteppingAction.o : $(LSExpSteppingAction_cc_dependencies)
	$(cpp_echo) $(src)LSExpSteppingAction.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpSteppingAction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpSteppingAction_cppflags) $(LSExpSteppingAction_cc_cppflags)  $(src)LSExpSteppingAction.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)LSExpDetectorConstruction.d

$(bin)$(binobj)LSExpDetectorConstruction.d :

$(bin)$(binobj)LSExpDetectorConstruction.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)LSExpDetectorConstruction.o : $(src)LSExpDetectorConstruction.cc
	$(cpp_echo) $(src)LSExpDetectorConstruction.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpDetectorConstruction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpDetectorConstruction_cppflags) $(LSExpDetectorConstruction_cc_cppflags)  $(src)LSExpDetectorConstruction.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(LSExpDetectorConstruction_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)LSExpDetectorConstruction.cc

$(bin)$(binobj)LSExpDetectorConstruction.o : $(LSExpDetectorConstruction_cc_dependencies)
	$(cpp_echo) $(src)LSExpDetectorConstruction.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpDetectorConstruction_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpDetectorConstruction_cppflags) $(LSExpDetectorConstruction_cc_cppflags)  $(src)LSExpDetectorConstruction.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)LSExpParticleGun.d

$(bin)$(binobj)LSExpParticleGun.d :

$(bin)$(binobj)LSExpParticleGun.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)LSExpParticleGun.o : $(src)LSExpParticleGun.cc
	$(cpp_echo) $(src)LSExpParticleGun.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpParticleGun_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpParticleGun_cppflags) $(LSExpParticleGun_cc_cppflags)  $(src)LSExpParticleGun.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(LSExpParticleGun_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)LSExpParticleGun.cc

$(bin)$(binobj)LSExpParticleGun.o : $(LSExpParticleGun_cc_dependencies)
	$(cpp_echo) $(src)LSExpParticleGun.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpParticleGun_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpParticleGun_cppflags) $(LSExpParticleGun_cc_cppflags)  $(src)LSExpParticleGun.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)LSExpPhysicsList.d

$(bin)$(binobj)LSExpPhysicsList.d :

$(bin)$(binobj)LSExpPhysicsList.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)LSExpPhysicsList.o : $(src)LSExpPhysicsList.cc
	$(cpp_echo) $(src)LSExpPhysicsList.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpPhysicsList_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpPhysicsList_cppflags) $(LSExpPhysicsList_cc_cppflags)  $(src)LSExpPhysicsList.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(LSExpPhysicsList_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)LSExpPhysicsList.cc

$(bin)$(binobj)LSExpPhysicsList.o : $(LSExpPhysicsList_cc_dependencies)
	$(cpp_echo) $(src)LSExpPhysicsList.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(LSExpPhysicsList_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(LSExpPhysicsList_cppflags) $(LSExpPhysicsList_cc_cppflags)  $(src)LSExpPhysicsList.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)DetPrototypeSvc.d

$(bin)$(binobj)DetPrototypeSvc.d :

$(bin)$(binobj)DetPrototypeSvc.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)DetPrototypeSvc.o : $(src)DetPrototypeSvc.cc
	$(cpp_echo) $(src)DetPrototypeSvc.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(DetPrototypeSvc_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(DetPrototypeSvc_cppflags) $(DetPrototypeSvc_cc_cppflags)  $(src)DetPrototypeSvc.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(DetPrototypeSvc_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)DetPrototypeSvc.cc

$(bin)$(binobj)DetPrototypeSvc.o : $(DetPrototypeSvc_cc_dependencies)
	$(cpp_echo) $(src)DetPrototypeSvc.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(DetPrototypeSvc_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(DetPrototypeSvc_cppflags) $(DetPrototypeSvc_cc_cppflags)  $(src)DetPrototypeSvc.cc

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (,)

ifneq ($(MAKECMDGOALS),DetSimOptionsclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)GlobalGeomInfo.d

$(bin)$(binobj)GlobalGeomInfo.d :

$(bin)$(binobj)GlobalGeomInfo.o : $(cmt_final_setup_DetSimOptions)

$(bin)$(binobj)GlobalGeomInfo.o : $(src)GlobalGeomInfo.cc
	$(cpp_echo) $(src)GlobalGeomInfo.cc
	$(cpp_silent) $(cppcomp)  -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(GlobalGeomInfo_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(GlobalGeomInfo_cppflags) $(GlobalGeomInfo_cc_cppflags)  $(src)GlobalGeomInfo.cc
endif
endif

else
$(bin)DetSimOptions_dependencies.make : $(GlobalGeomInfo_cc_dependencies)

$(bin)DetSimOptions_dependencies.make : $(src)GlobalGeomInfo.cc

$(bin)$(binobj)GlobalGeomInfo.o : $(GlobalGeomInfo_cc_dependencies)
	$(cpp_echo) $(src)GlobalGeomInfo.cc
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DetSimOptions_pp_cppflags) $(lib_DetSimOptions_pp_cppflags) $(GlobalGeomInfo_pp_cppflags) $(use_cppflags) $(DetSimOptions_cppflags) $(lib_DetSimOptions_cppflags) $(GlobalGeomInfo_cppflags) $(GlobalGeomInfo_cc_cppflags)  $(src)GlobalGeomInfo.cc

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: DetSimOptionsclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(DetSimOptions.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

DetSimOptionsclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library DetSimOptions
	-$(cleanup_silent) cd $(bin) && \rm -f $(library_prefix)DetSimOptions$(library_suffix).a $(library_prefix)DetSimOptions$(library_suffix).$(shlibsuffix) DetSimOptions.stamp DetSimOptions.shstamp
#-- end of cleanup_library ---------------
