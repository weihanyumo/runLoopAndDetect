//
//  ViewController.m
//  detectDemo
//
//  Created by renxin on 15/6/10.
//  Copyright (c) 2015年 dhd. All rights reserved.
//

#import "ViewController.h"
#import <opencv2/opencv.hpp>
#import <opencv2/imgproc/types_c.h>
#import <opencv2/imgcodecs/ios.h>
#import <opencv2/features2d.hpp>
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#import "UpixDBSCAN.h"
#import <CoreAudio/CoreAudioTypes.h>
#import <opencv2/core/types_c.h>
#import <dispatch/dispatch.h>
#include <CoreFoundation/CoreFoundation.h>
#include <sys/socket.h>
#include <netinet/in.h>
#import <mach/mach_time.h>
#import <CoreFoundation/CFRunLoop.h>



BOOL threadProcess1Finished =NO;
BOOL threadProcess2Finished =NO;
BOOL pageStillLoading = YES;

@interface ViewController ()
{
    IBOutlet UIButton        *btnNormal;
    IBOutlet UIButton        *btnRunloop;
    IBOutlet UIButton        *btnTest;
    IBOutlet UIButton        *btnObserver;
    IBOutlet UIButton        *btnTestLeak;
}
@property(nonatomic, strong) NSMutableArray *arrayLeak;
@end

@implementation ViewController

- (void)viewDidLoad {

    [super viewDidLoad];
    
//    createCustomSource();
    //createPortSource();
//    [self testOpenCV];
}

//
-(IBAction)btnTestLeakClicked:(id)sender
{
    self.arrayLeak = [[NSMutableArray alloc]initWithCapacity:0];
    for(int i=0;i<1000;i++)
    {
//        int *leak = (int*)malloc(sizeof(int) * 100);
//        NSLog([NSString stringWithFormat:@"%d", leak[0]]);
        [_arrayLeak addObject:@"dadafd;"];
    }
}

//opencv
-(void) testOpenCV{
    cv::Mat cvImage;
    UIImage *image = [UIImage imageNamed:@"1.png"];
    UIImage *image2 = [UIImage imageNamed:@"2.png"];
    // Convert UIImage * to cv::Mat
    //UIImageToMat(image, cvImage);
    [self matchImage:image And:image2];
    
    return;
    if (!cvImage.empty())
    {
        cv::Mat gray;
        cv::cvtColor(cvImage, gray, CV_RGBA2GRAY);
        cv::GaussianBlur(gray, gray, cv::Size(5,5), 1.2,1.2);
        cv::Mat edges;
        cv::Canny(gray, edges, 0, 60);
        cvImage.setTo(cv::Scalar::all(255));
        cvImage.setTo(cv::Scalar(0,128,255,255),edges);
        
        UIImageView *imgView = [[UIImageView alloc] initWithImage:MatToUIImage(cvImage)];
        
        [self.view addSubview:imgView];
    }
}
- (BOOL) matchImage:(UIImage *)image1 And:(UIImage*)image2{
    cv::Mat img1;
    cv::Mat img2;
    
    UIImageToMat(image1, img1, YES);
    UIImageToMat(image2, img2, YES);
    
    std::vector<cv::KeyPoint> key_points_1;
    std::vector<cv::KeyPoint> key_points_2;
    cv::Mat descriptors1;
    cv::Mat descriptors2;
    std::vector<cv::DMatch> matches;
    
    cv::Ptr<cv::FastFeatureDetector> detector = cv::FastFeatureDetector::create();
    cv::ORB *orb = cv::ORB::create(400, 1.2, 8, 21, 0, 2, cv::ORB::FAST_SCORE, 21,20);
//    cv::ORB *orb2 = cv::ORB::create(400, 1.2, 8, 21, 0, 2, cv::ORB::FAST_SCORE, 21,20);
//    cv::Ptr<cv::DescriptorExtractor>descriptorMatcher2 = cv::DescriptorExtractor::create("");
    cv::Ptr<cv::DescriptorMatcher> descriptorMatcher = cv::DescriptorMatcher::create("BruteForce-Hamming");
#if 0
    orb->detectAndCompute(img1, cv::Mat(), key_points_1, descriptors1);
    orb->detectAndCompute(img2, cv::Mat(), key_points_2, descriptors2);
    
    orb->detect(img1, key_points_1);
    orb2->detect(img2, key_points_2);
    orb->compute(img1, key_points_1, descriptors1);
    orb2->compute(img2, key_points_2, descriptors2);
#else
    detector->detect(img1, key_points_1);
    detector->detect(img2, key_points_2);
//    orb->compute(img1, key_points_1, descriptors1);
//    orb->compute(img2, key_points_2, descriptors2);
    
    detector->compute(img1, key_points_1, descriptors1);
    detector->compute(img2, key_points_2, descriptors2);
#endif
    descriptorMatcher->match( descriptors1, descriptors2, matches);
    
    
    double max_dist = 0; double min_dist = 100;
    for( int i = 0; i < descriptors1.rows; i++ )
    {
        double dist = matches[i].distance;
        if( dist < min_dist ) min_dist = dist;
        if( dist > max_dist ) max_dist = dist;
    }
    
    std::vector< cv::DMatch > good_matches;
    for( int i = 0; i < descriptors1.rows; i++ )
    {
        if( matches[i].distance < 0.6*max_dist )
        {
            good_matches.push_back( matches[i]);
        }
    }
    
    
    for (int i=0; i<key_points_1.size(); i++)
    {
        cv::KeyPoint keypoint = key_points_1.at(i);
        cv::Point ipt = keypoint.pt;
        cv::circle(img1, ipt, 2, cv::Scalar(255, 0,0));

    }
    
    
    UIImageView *imgView = [[UIImageView alloc] initWithImage:MatToUIImage(img1)];
    [self.view addSubview:imgView];

    return NO;
}

//GCD
- (void) testGroup{
    NSLog(@"....................testGroup...........................");
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
    dispatch_queue_t myQueue = dispatch_queue_create("test.myQueue", nil);
    dispatch_group_t group = dispatch_group_create();
    for(int i=0; i<10; i++)
    {
        dispatch_group_async(group, myQueue/*queue*/, ^{
            NSLog(@"gorup %d", i);
        });
    }
    dispatch_group_notify(group, queue, ^{
        NSLog(@"group end!!!");
    });
}
- (void) testApply{
    NSLog(@"....................testApply...........................");
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
    dispatch_apply(10, queue, ^(size_t index){
       NSLog(@"apply %zu", index);
    });
}
- (void) testSource{
    NSLog(@"....................testSource...........................");
    dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
    dispatch_queue_t myQueue = dispatch_queue_create("test.test.myQueue", nil);
    dispatch_source_t source = dispatch_source_create(DISPATCH_SOURCE_TYPE_DATA_ADD, 0, 0, dispatch_get_main_queue());
    dispatch_source_set_event_handler(source, ^{
        NSLog(@"current progress:%ld", dispatch_source_get_data(source));
    });
    dispatch_resume(source);
    
    dispatch_apply(20, globalQueue, ^(size_t index) {
        NSLog(@"apply %ld", index);
        dispatch_source_merge_data(source, 1);
        
    });
}


//test runloop
-(void)threadProce1{
    NSLog(@"Enter threadProce1.");
    for (int i=0; i<5;i++)
    {
        NSLog(@"InthreadProce1 count = %d.", i);
        sleep(1);
    }
    threadProcess1Finished =YES;
    NSLog(@"Exit threadProce1.");
}
-(void)threadProce2{
    NSLog(@"Enter threadProce2.");
    for (int i=0; i<5;i++)
    {
        NSLog(@"InthreadProce2 count = %d.", i);
        sleep(1);
    }
    threadProcess2Finished =YES;
    NSLog(@"Exit threadProce2.");
}
-(void)threadProce3{
    [self addObserverToCurrentRunloop];
    
    NSLog(@"Enter threadProce3.");
    for (int i=0; i<5;i++)
    {
        NSLog(@"InthreadProce3 count = %d.", i);
        sleep(1);
    }
    NSLog(@"Exit threadProce3.");
    [self createCustomSource];
}

- (IBAction)buttonNormalThreadTestPressed:(UIButton *)sender {
    NSLog(@"EnterbuttonNormalThreadTestPressed");
    threadProcess1Finished =NO;
    NSLog(@"Start a new thread.");
    [NSThread detachNewThreadSelector: @selector(threadProce1)
                            toTarget: self
                           withObject: nil];
    
    while (!threadProcess1Finished)
    {
        [NSThread sleepForTimeInterval: 0.5];
    }
    NSLog(@"ExitbuttonNormalThreadTestPressed");
}
- (IBAction)buttonRunloopPressed:(id)sender {
    NSLog(@"Enter buttonRunloopPressed");
    threadProcess2Finished =NO;
    NSLog(@"Start a new thread.");
    [NSThread detachNewThreadSelector: @selector(threadProce2)
                            toTarget: self
                          withObject: nil];
    
    while (!threadProcess2Finished)
    {
//        NSLog(@"Begin runloop");
        [[NSRunLoop currentRunLoop] runMode:NSDefaultRunLoopMode beforeDate: [NSDate distantFuture]];
//        NSLog(@"End runloop.");
    }
    NSLog(@"Exit buttonRunloopPressed");
}
- (IBAction)buttonTestPressed:(id)sender{
    NSLog(@"Enter buttonTestPressed");
    NSLog(@"Exit buttonTestPressed");
    
    
}
- (IBAction)buttonTestObserver:(id)sender{
    [NSThread detachNewThreadSelector: @selector(threadProce3)
                             toTarget: self
                           withObject: nil];
}

//
void myRunLoopObserver(CFRunLoopObserverRef observer, CFRunLoopActivity activity, void *info)
{
    NSLog(@"activity is %lu", activity);
}

CFDataRef myCallbackFunc(CFMessagePortRef local, SInt32 msgid, CFDataRef data, void *info)
{
    return nil;
}

-(void) createCustomSource
{
    CFRunLoopSourceContext context = {0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    CFRunLoopSourceRef source = CFRunLoopSourceCreate(kCFAllocatorDefault, 0, &context);
    CFRunLoopAddSource(CFRunLoopGetCurrent(), source, kCFRunLoopDefaultMode);
    while (pageStillLoading)
    {
        CFRunLoopRun();
    }
    CFRunLoopRemoveSource(CFRunLoopGetCurrent(), source, kCFRunLoopDefaultMode);
    CFRelease(source);
}
void createPortSource()
{
    CFMessagePortRef port = CFMessagePortCreateLocal(kCFAllocatorDefault, CFSTR("com.someport"),myCallbackFunc, NULL, NULL);
    CFRunLoopSourceRef source =  CFMessagePortCreateRunLoopSource(kCFAllocatorDefault, port, 0);
    CFRunLoopAddSource(CFRunLoopGetCurrent(), source, kCFRunLoopCommonModes);
    while (pageStillLoading)
    {
        CFRunLoopRun();
    }
    CFRunLoopRemoveSource(CFRunLoopGetCurrent(), source, kCFRunLoopDefaultMode);
    CFRelease(source);
}

- (void) addObserverToCurrentRunloop
{
    NSRunLoop*myRunLoop = [NSRunLoop currentRunLoop];
    CFRunLoopObserverContext  context = {0, (__bridge void*)self, NULL, NULL, NULL};
    CFRunLoopObserverRef    observer = CFRunLoopObserverCreate(kCFAllocatorDefault,
                                                              kCFRunLoopBeforeTimers, YES, 0, &myRunLoopObserver, &context);
    if (observer)
    {
        CFRunLoopRef    cfLoop = [myRunLoop getCFRunLoop];
       // cfLoop = CFRunLoopGetMain();
        CFRunLoopAddObserver(cfLoop, observer, kCFRunLoopDefaultMode);
    }
}

@end
