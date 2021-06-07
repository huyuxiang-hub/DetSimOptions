# echo "setup DetSimOptions v0 in /afs/ihep.ac.cn/users/h/huyuxiang/junofs/juno_centos7/offline/Simulation/DetSimV2"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J21v1r0-branch-python3/ExternalLibs/CMT/v1r26
endif
source ${CMTROOT}/mgr/setup.csh
set cmtDetSimOptionstempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtDetSimOptionstempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=DetSimOptions -version=v0 -path=/afs/ihep.ac.cn/users/h/huyuxiang/junofs/juno_centos7/offline/Simulation/DetSimV2  -no_cleanup $* >${cmtDetSimOptionstempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=DetSimOptions -version=v0 -path=/afs/ihep.ac.cn/users/h/huyuxiang/junofs/juno_centos7/offline/Simulation/DetSimV2  -no_cleanup $* >${cmtDetSimOptionstempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtDetSimOptionstempfile}
  unset cmtDetSimOptionstempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtDetSimOptionstempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtDetSimOptionstempfile}
unset cmtDetSimOptionstempfile
exit $cmtsetupstatus

