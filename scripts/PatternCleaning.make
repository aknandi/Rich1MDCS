#-- start of make_header -----------------

#====================================
#  Application ../scripts/PatternCleaning
#
#   Generated Fri Apr  1 17:13:22 2016  by anandi
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_../scripts/PatternCleaning_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_../scripts/PatternCleaning_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_../scripts/PatternCleaning

Rich1MDCS_tag = $(tag)

#cmt_local_tagfile_../scripts/PatternCleaning = $(Rich1MDCS_tag)_../scripts/PatternCleaning.make
cmt_local_tagfile_../scripts/PatternCleaning = $(bin)$(Rich1MDCS_tag)_../scripts/PatternCleaning.make

else

tags      = $(tag),$(CMTEXTRATAGS)

Rich1MDCS_tag = $(tag)

#cmt_local_tagfile_../scripts/PatternCleaning = $(Rich1MDCS_tag).make
cmt_local_tagfile_../scripts/PatternCleaning = $(bin)$(Rich1MDCS_tag).make

endif

include $(cmt_local_tagfile_../scripts/PatternCleaning)
#-include $(cmt_local_tagfile_../scripts/PatternCleaning)

ifdef cmt_../scripts/PatternCleaning_has_target_tag

cmt_final_setup_../scripts/PatternCleaning = $(bin)setup_../scripts/PatternCleaning.make
#cmt_final_setup_../scripts/PatternCleaning = $(bin)Rich1MDCS_../scripts/PatternCleaningsetup.make
cmt_local_../scripts/PatternCleaning_makefile = $(bin)../scripts/PatternCleaning.make

else

cmt_final_setup_../scripts/PatternCleaning = $(bin)setup.make
#cmt_final_setup_../scripts/PatternCleaning = $(bin)Rich1MDCSsetup.make
cmt_local_../scripts/PatternCleaning_makefile = $(bin)../scripts/PatternCleaning.make

endif

cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)Rich1MDCSsetup.make

#../scripts/PatternCleaning :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) '../scripts/PatternCleaning'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = ../scripts/PatternCleaning/
#../scripts/PatternCleaning::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of application_header

../scripts/PatternCleaning :: dirs  $(bin)../scripts/PatternCleaning${application_suffix}
	$(echo) "../scripts/PatternCleaning ok"

#-- end of application_header
#-- start of application

$(bin)../scripts/PatternCleaning${application_suffix} :: $(bin)CleanRun.o $(use_stamps) $(../scripts/PatternCleaning_stamps) $(../scripts/PatternCleaningstamps) $(use_requirements)
	$(link_echo) "application $@"
	$(link_silent) $(cpplink) -o $(@).new $(bin)CleanRun.o $(cmt_installarea_linkopts) $(../scripts/PatternCleaning_use_linkopts) $(../scripts/PatternCleaninglinkopts) && mv -f $(@).new $(@)

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

install_dir = ${CMTINSTALLAREA}/$(tag)/bin
../scripts/PatternCleaninginstallname = ../scripts/PatternCleaning${application_suffix}

../scripts/PatternCleaning :: ../scripts/PatternCleaninginstall

install :: ../scripts/PatternCleaninginstall

../scripts/PatternCleaninginstall :: $(install_dir)/$(../scripts/PatternCleaninginstallname)
ifdef CMTINSTALLAREA
	$(echo) "installation done"
endif

$(install_dir)/$(../scripts/PatternCleaninginstallname) :: $(bin)$(../scripts/PatternCleaninginstallname)
ifdef CMTINSTALLAREA
	$(install_silent) $(cmt_install_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(../scripts/PatternCleaninginstallname)" \
	    -out "$(install_dir)" \
	    -cmd "$(cmt_installarea_command)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

##../scripts/PatternCleaningclean :: ../scripts/PatternCleaninguninstall

uninstall :: ../scripts/PatternCleaninguninstall

../scripts/PatternCleaninguninstall ::
ifdef CMTINSTALLAREA
	$(cleanup_silent) $(cmt_uninstall_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(../scripts/PatternCleaninginstallname)" \
	    -out "$(install_dir)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

#	@echo "------> (../scripts/PatternCleaning.make) Removing installed files"
#-- end of application
#-- start of dependency ------------------
ifneq ($(MAKECMDGOALS),../scripts/PatternCleaningclean)
ifneq ($(MAKECMDGOALS),uninstall)

#$(bin)../scripts/PatternCleaning_dependencies.make :: dirs

ifndef QUICK
$(bin)../scripts/PatternCleaning_dependencies.make : ../scripts/CleanRun.C $(use_requirements) $(cmt_final_setup_../scripts/PatternCleaning)
	$(echo) "(../scripts/PatternCleaning.make) Rebuilding $@"; \
	  $(build_dependencies) ../scripts/PatternCleaning -all_sources -out=$@ ../scripts/CleanRun.C
endif

#$(../scripts/PatternCleaning_dependencies)

-include $(bin)../scripts/PatternCleaning_dependencies.make

endif
endif
#-- end of dependency -------------------
#-- start of cpp ------

ifneq ($(MAKECMDGOALS),../scripts/PatternCleaningclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)CleanRun.d
endif
endif


$(bin)$(binobj)CleanRun.o $(bin)$(binobj)CleanRun.d : ../scripts/CleanRun.C  $(use_requirements) $(cmt_final_setup_../scripts/PatternCleaning)
	$(cpp_echo) ../scripts/CleanRun.C
	@mkdir -p $(@D)
	$(cpp_silent) $(cppcomp) $(use_pp_cppflags) $(../scripts/PatternCleaning_pp_cppflags) $(app_../scripts/PatternCleaning_pp_cppflags) $(CleanRun_pp_cppflags) $(use_cppflags) $(../scripts/PatternCleaning_cppflags) $(lib_../scripts/PatternCleaning_cppflags) $(app_../scripts/PatternCleaning_cppflags) $(CleanRun_cppflags) $(CleanRun_C_cppflags) -I../scripts -MP -MMD -MT $(bin)$(binobj)CleanRun.o -MT $(bin)$(binobj)CleanRun.d -MF $(bin)$(binobj)CleanRun.d -o $(bin)$(binobj)CleanRun.o ../scripts/CleanRun.C


#-- end of cpp ------
#-- start of cleanup_header --------------

clean :: ../scripts/PatternCleaningclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(../scripts/PatternCleaning.make) $@: No rule for such target" >&2
#	@echo "#CMT> Warning: $@: No rule for such target" >&2; exit
	if echo $@ | grep '$(package)setup\.make$$' >/dev/null; then\
	 echo "$(CMTMSGPREFIX)" "(../scripts/PatternCleaning.make): $@: File no longer generated" >&2; exit 0; fi
else
.DEFAULT::
	$(echo) "(../scripts/PatternCleaning.make) PEDANTIC: $@: No rule for such target" >&2
	if echo $@ | grep '$(package)setup\.make$$' >/dev/null; then\
	 echo "$(CMTMSGPREFIX)" "(../scripts/PatternCleaning.make): $@: File no longer generated" >&2; exit 0;\
	 elif test $@ = "$(cmt_final_setup)" -o\
	 $@ = "$(cmt_final_setup_../scripts/PatternCleaning)" ; then\
	 found=n; for s in 1 2 3 4 5; do\
	 sleep $$s; test ! -f $@ || { found=y; break; }\
	 done; if test $$found = n; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(../scripts/PatternCleaning.make) PEDANTIC: $@: Seems to be missing. Ignore it for now" >&2; exit 0 ; fi;\
	 elif test `expr $@ : '.*/'` -ne 0 ; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(../scripts/PatternCleaning.make) PEDANTIC: $@: Seems to be a missing file. Please check" >&2; exit 2 ; \
	 else\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(../scripts/PatternCleaning.make) PEDANTIC: $@: Seems to be a fake target due to some pattern. Just ignore it" >&2 ; exit 0; fi
endif

../scripts/PatternCleaningclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_application ------
	$(cleanup_echo) application ../scripts/PatternCleaning
	-$(cleanup_silent) cd $(bin); /bin/rm -f ../scripts/PatternCleaning${application_suffix}
#-- end of cleanup_application ------
#-- start of cleanup_objects ------
	$(cleanup_echo) objects ../scripts/PatternCleaning
	-$(cleanup_silent) /bin/rm -f $(bin)CleanRun.o
	-$(cleanup_silent) /bin/rm -f $(patsubst %.o,%.d,$(bin)CleanRun.o) $(patsubst %.o,%.dep,$(bin)CleanRun.o) $(patsubst %.o,%.d.stamp,$(bin)CleanRun.o)
	-$(cleanup_silent) cd $(bin); /bin/rm -rf ../scripts/PatternCleaning_deps ../scripts/PatternCleaning_dependencies.make
#-- end of cleanup_objects ------
