import cv2
import numpy as np
import pandas as pd
from os import path
from os import mkdir
from sys import argv
import time
#script, basedir, start_frame_str, end_frame_str, name = argv

if __name__ == "__main__":
    # TODO Tim: Replace this with command line arguments and argparse when sync works.
    #output_path = path.join(r'SyncedWithTimestamps')
    output_path = 'D:\Log_D2P282649_2019.05.30_16.46.20\SyncedWithTimestamps_1'
    print('output path {} exists {}'.format(output_path,path.isdir(output_path)))
    #video_0_path = path.join(r'EndoscopeImageMemory_0.avi')
    #video_0_path = path.join(basedir,r'EndoscopeImageMemory_0_deint.avi')
    video_0_path = 'D:\Log_D2P282649_2019.05.30_16.46.20\part0001\EndoscopeImageMemory_0_deint.avi'
    #video_1_path = path.join(r'EndoscopeImageMemory_1.avi')
    #video_1_path = path.join(basedir,r'EndoscopeImageMemory_1_deint.avi')
    video_1_path = 'D:\Log_D2P282649_2019.05.30_16.46.20\part0001\EndoscopeImageMemory_1_deint.avi'
    #video_csv_0_path = path.join(r'EndoscopeImageMemory_0.csv')
    #video_csv_0_path = path.join(basedir,r'EndoscopeImageMemory_0.deint.csv')
    video_csv_0_path = 'D:\Log_D2P282649_2019.05.30_16.46.20\part0001\EndoscopeImageMemory_0.csv'
    #video_csv_1_path = path.join(r'EndoscopeImageMemory_1.csv')
    #video_csv_1_path = path.join(basedir,r'EndoscopeImageMemory_1.deint.csv')
    video_csv_1_path = 'D:\Log_D2P282649_2019.05.30_16.46.20\part0001\EndoscopeImageMemory_1.csv'
    #kinematics_csv_path = path.join(r'DaVinciSiMemory.csv')
    #kinematics_csv_path = 'D:\FFOutput\DaVinciSiMemory.csv'

    video_2_path = 'D:\Log_D2P282649_2019.05.30_16.46.20\part0001\RenderedImageMemory_0.avi'
    video_csv_2_path = 'D:\Log_D2P282649_2019.05.30_16.46.20\part0001\RenderedImageMemory_0.csv'

    if path.isdir(output_path) == False:
        mkdir(output_path)
    # video_0_path = 'E:\C:\Data\Images\TimK\ColonPhantomWith_daVinci\part0002\EndoscopeImageMemory_0.avi'
    # video_1_path = 'E:\\Software\\UCL\\daVinci\\MatLabScripts\\Data\\KnotTying\\Expert\\Log_D2P280782_2017.08.30_16.39.01_part0012\\EndoscopeImageMemory_1_768_576.avi'
	
    # Extract `header.timestamp` field for each video CSV
    timestamps_0 = pd.read_csv(video_csv_0_path)['header.timestamp']
    timestamps_1 = pd.read_csv(video_csv_1_path)['header.timestamp']
    timestamps_2 = pd.read_csv(video_csv_2_path)['header.timestamp']

    # Load the whole kinematics CSV and store `header.timestamp` separately
    #kinematics = pd.read_csv(kinematics_csv_path)
    #kinematics_timestamps = kinematics['header.timestamp']
    start_frame = int(1)
    end_frame = int(85680)

    #
    framerate = 30.0
    mus_per_frame = 1000.0 / framerate
    if start_frame < 0 :
        start_frame = 0
    # start_time = max(timestamps_0[0], timestamps_1[0])
    start_time = timestamps_2[start_frame]
    # end_time = min(timestamps_0[len(timestamps_0) - 1], timestamps_1[len(timestamps_1) - 1])
    if end_frame >= timestamps_2.size:
        end_frame = timestamps_2.size-1

    end_time = timestamps_2[end_frame]
    # Prepare containers for frame indices
    total_frames = int(np.ceil((end_time - start_time) / mus_per_frame))+1
    frames_0 = np.zeros((total_frames,))
    frames_1 = np.zeros((total_frames,))
    frames_2 = np.zeros((total_frames,))    
    #kinematics_indices = np.zeros((total_frames,))
    print('Final videos will have {} frames.'.format(total_frames))

    #
    frame_index = 0
    index_0 = int(0)
    index_1 = int(0)
    index_2 = int(0)
    #kinematics_index = int(0)
    #kinematics_offset = -130
    current_time = start_time
    print('Synchronizing frame indices...')
    while current_time < end_time:

        # Find closest frame for first video
        while True:
            last_frame_diff = np.abs(current_time - timestamps_0[index_0])
            next_index = index_0 + 1
            if next_index >= len(timestamps_0):
                break
            frame_diff = np.abs(current_time - timestamps_0[next_index])
            if frame_diff <= last_frame_diff:
                index_0 = next_index
            else:
                break

        # Find closest frame for second video
        while True:
            last_frame_diff = np.abs(current_time - timestamps_1[index_1])
            next_index = index_1 + 1
            if next_index >= len(timestamps_1):
                break
            frame_diff = np.abs(current_time - timestamps_1[next_index])
            if frame_diff <= last_frame_diff:
                index_1 = next_index
            else:
                break

        # Find closest frame for third video
        while True:
            last_frame_diff = np.abs(current_time - timestamps_2[index_2])
            next_index = index_2 + 1
            if next_index >= len(timestamps_2):
                break
            frame_diff = np.abs(current_time - timestamps_2[next_index])
            if frame_diff <= last_frame_diff:
                index_2 = next_index
            else:
                break

        # # # Find closest kinematics frame
        #while True:
        #    last_frame_diff = np.abs(current_time - kinematics_timestamps[kinematics_index]+kinematics_offset)
        #    next_index = kinematics_index + 1
        #    if next_index >= len(kinematics_timestamps):
        #        break
        #    frame_diff = np.abs(current_time - kinematics_timestamps[next_index]+kinematics_offset)
        #    if frame_diff < last_frame_diff:
        #        kinematics_index = next_index
        #    else:
        #        break

        # Store corresponding frames
        frames_0[frame_index] = index_0
        frames_1[frame_index] = index_1
        frames_2[frame_index] = index_2
        #kinematics_indices[frame_index] = kinematics_index
        current_time += mus_per_frame
        frame_index += 1
        print('Frame_idx={}, end_time-current_time={}'.format(frame_index,end_time-current_time))
    print('Frame_idx={}, total_frames={}'.format(frame_index,total_frames))
    if frame_index == total_frames-1 :
        total_frames = frame_index
    assert frame_index == total_frames
    print('Done!')
    print('frames_0[0], frames_1[0], frames_2[0] = {},{},{}'.format(frames_0[0], frames_1[0], frames_2[0]))
    print('frames_0[total_frames-1], frames_1[total_frames-1], frames_2[total_frames-1] = {},{},{}'.format(
       frames_0[total_frames-1], frames_1[total_frames-1], frames_2[total_frames-1]))
    print('Timestamp difference frames_0[end].timestamp-frames_0[start].timestamp = {}'.format(
       timestamps_0[frames_0[total_frames-1]] - timestamps_0[frames_0[0]]))
    print('Timestamp difference frames_1[end].timestamp-frames_1[start].timestamp = {}'.format(
       timestamps_1[frames_1[total_frames-1]] - timestamps_1[frames_1[0]]))
    print('Timestamp difference frames_2[end].timestamp-frames_2[start].timestamp = {}'.format(
       timestamps_2[frames_2[total_frames-1]] - timestamps_2[frames_2[0]]))       
    #print('Timestamp difference kinematics[end].timestamp-kinematics[start].timestamp = {}'.format(
       #kinematics_timestamps[kinematics_indices[total_frames - 1]] - kinematics_timestamps[kinematics_indices[0]]))
    print('frames_0[100], frames_1[100], frames_2[100] = {},{},{}'.format(frames_0[100], frames_1[100], frames_2[100]))
    # Write out the synced kinematics CSV
    #kinematics_output_path = path.join(output_path, path.basename(kinematics_csv_path))
    #base = path.splitext(kinematics_output_path)[0]
    #kinematics_output_path = base + '_test.csv'
    #kinematics_output_dataframe = pd.DataFrame(columns=kinematics.columns.values)
    #print('Kinematics CSV will be written to `{}`...'.format(kinematics_output_path))
    #print('Preparing kinematics CSV...')
    #for i in range(total_frames):
    #    print('Processing row {} out of {}...'.format(i + 1, total_frames), end='\r')
    #    kin_index = kinematics_indices[i]
    #    kinematics_output_dataframe.loc[i] = kinematics.loc[kin_index]
    #kinematics_output_dataframe.to_csv(kinematics_output_path, header=True)
    print('Done!')

    video_0_source = cv2.VideoCapture(video_0_path)
    if (video_0_source.isOpened() == False): 
        print("Error opening video file 0")
    video_1_source = cv2.VideoCapture(video_1_path)
    if (video_1_source.isOpened() == False): 
        print("Error opening video file 1")
    video_2_source = cv2.VideoCapture(video_2_path)
    if (video_2_source.isOpened() == False): 
        print("Error opening video file 2")

    width = video_0_source.get(cv2.CAP_PROP_FRAME_WIDTH)	
    height = video_0_source.get(cv2.CAP_PROP_FRAME_HEIGHT)	
    width2 = video_2_source.get(cv2.CAP_PROP_FRAME_WIDTH)	
    height2 = video_2_source.get(cv2.CAP_PROP_FRAME_HEIGHT)	
    video_0_output_path = path.join(output_path, path.basename(video_0_path))
    base = path.splitext(video_0_output_path)[0]
    video_0_output_path = base + '_test.avi'
    video_1_output_path = path.join(output_path, path.basename(video_1_path))
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    base = path.splitext(video_1_output_path)[0]
    video_1_output_path = base + '_test.avi'
    video_2_output_path = path.join(output_path, path.basename(video_2_path))
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    base = path.splitext(video_2_output_path)[0]
    video_2_output_path = base + '_test.avi'

    video_0_output = cv2.VideoWriter(video_0_output_path, fourcc, framerate, (int(width), int(height)))
    video_1_output = cv2.VideoWriter(video_1_output_path, fourcc, framerate, (int(width), int(height)))
    video_2_output = cv2.VideoWriter(video_2_output_path, fourcc, framerate, (int(width2), int(height2)))

    print('Videos will be written to `{}`, `{}` and `{}`...'.format(video_0_output_path, video_1_output_path, video_2_output_path))

    print('Preparing videos...')
    ret, last_frame_0 = video_0_source.read()
    print('Video 0 Framesize is {} by {}',video_0_source.get(3), video_0_source.get(4))
    #gray = cv2.cvtColor(last_frame_0, cv2.COLOR_BGR2GRAY)
    #cv2.imshow('last_frame_0',last_frame_0)
    ret, last_frame_1 = video_1_source.read()
    #gray = cv2.cvtColor(last_frame_1, cv2.COLOR_BGR2GRAY)
    #cv2.imshow('last_frame_1',last_frame_1)
    ret, last_frame_2 = video_2_source.read()

    last_frame_index_0 = 0
    last_frame_index_1 = 0
    last_frame_index_2 = 0
    for i in range(total_frames):

        target_frame_index_0 = frames_0[i]
        target_frame_index_1 = frames_1[i]
        target_frame_index_2 = frames_2[i]

        print('Processing frame {} out of {} - target frames are {}, {} and {}'.format(i + 1, total_frames, target_frame_index_0, target_frame_index_1, target_frame_index_2), end='\n')

        # Skip until the desired frame in the first video
        if last_frame_index_0 != target_frame_index_0:
            backup_frame = None
            while last_frame_index_0 != target_frame_index_0 - 1:
                last_frame_index_0 += 1
                success, frame = video_0_source.read()
                if success:
                    backup_frame = frame

            last_frame_index_0 += 1
            success, frame = video_0_source.read()
            if success:
                last_frame_0 = frame
            elif backup_frame is not None:
                last_frame_0 = backup_frame

        # Skip until the desired frame in the second video
        if last_frame_index_1 != target_frame_index_1:
            backup_frame = None
            while last_frame_index_1 != target_frame_index_1 - 1:
                last_frame_index_1 += 1
                success, frame = video_1_source.read()
                if success:
                    backup_frame = frame

            last_frame_index_1 += 1
            success, frame = video_1_source.read()
            if success:
                last_frame_1 = frame
            elif backup_frame is not None:
                last_frame_1 = backup_frame

            # Skip until the desired frame in the second video
        if last_frame_index_2 != target_frame_index_2:
            backup_frame = None
            while last_frame_index_2 != target_frame_index_2 - 1:
                last_frame_index_2 += 1
                success, frame = video_2_source.read()
                if success:
                    backup_frame = frame

            last_frame_index_2 += 1
            success, frame = video_2_source.read()
            if success:
                last_frame_2 = frame
            elif backup_frame is not None:
                last_frame_2 = backup_frame

        # Write out the current frame
        font  = cv2.FONT_HERSHEY_SIMPLEX
        bottomLeftCornerOfText = (int(width-450),35)
        fontScale              = 1
        fontColor              = (255,255,255)
        lineType               = 2
        readable = time.ctime(int(timestamps_0[last_frame_index_0])/1000)
        print('Writing video 0 frame {} out of {}...'.format(i + 1, total_frames), end='\n')
        # Deinterlace_Easy - doesn't work
        #last_frame_0[1::2]=last_frame_0[::2]
        #last_frame_0[1:-1:2] = (last_frame_0[1:-2:2]+last_frame_0[2::2])/2
        cv2.putText(last_frame_0,readable, bottomLeftCornerOfText, font, fontScale, fontColor, lineType)
        video_0_output.write(last_frame_0)
        readable = time.ctime(int(timestamps_1[last_frame_index_1])/1000)
        print('Writing video 1 frame {} out of {}...'.format(i + 1, total_frames), end='\n')
        # Deinterlace_Easy - doesn't work
        #last_frame_1[1::2]=last_frame_1[::2]
        #last_frame_1[1:-1:2] = (last_frame_1[1:-2:2]+last_frame_1[2::2])/2
        cv2.putText(last_frame_1,readable, bottomLeftCornerOfText, font, fontScale, fontColor, lineType)
        video_1_output.write(last_frame_1)
        print('Writing video 2 frame {} out of {}...'.format(i + 1, total_frames), end='\n')
        cv2.putText(last_frame_2,readable, bottomLeftCornerOfText, font, fontScale, fontColor, lineType)
        video_2_output.write(last_frame_2)
        print('i={}'.format(i))

    video_0_source.release()
    video_1_source.release()
    video_2_source.release()
    video_0_output.release()
    video_1_output.release()
    video_2_output.release()

    print('Done!')
    cv2.destroyAllWindows()
