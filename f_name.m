 % Creates file name. 1 = faces; 2 = places.

function [file_name1,file_name2] = f_name(ID,subject_id,run_name);

        %  file name example 'sub-01_run-1_timeseries.mat'
          if ID<10;
            sub_name = sprintf('%s%d','sub-0',subject_id(ID,1));
          else
            sub_name = sprintf('%s%d','sub-',subject_id(ID,1));
          end
    
        %run correspond to scanning session align to face onset times
          if run_name(ID,1)== 1;
              run_newname= 2;
          else run_name(ID,1)== 2;
              run_newname= 1;
          end
        
         run = sprintf('%s%d','run-',run_name(ID));
         run2 = sprintf('%s%d','run-',run_newname);
            
        %filename of time series
            file_name1 = sprintf('%s%s%s%s%s',sub_name,'_',run,'_','timeseries.mat');
            file_name2 = sprintf('%s%s%s%s%s',sub_name,'_',run2,'_','timeseries.mat');

        
            