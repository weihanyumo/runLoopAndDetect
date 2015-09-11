//
//  main.m
//  detectDemo
//
//  Created by renxin on 15/6/10.
//  Copyright (c) 2015年 dhd. All rights reserved.
//

#import "AppDelegate.h"

#if 1
int main(int argc, char * argv[]) {
    @autoreleasepool {
        return UIApplicationMain(argc, argv, nil, NSStringFromClass([AppDelegate class]));
    }
}
#else

static void _perform(void *info __unused)
{
    printf("hello\n");
}

static void _timer(CFRunLoopTimerRef timer __unused, void *info)
{
    CFRunLoopSourceSignal(info);
}

CFDataRef myCallbackFunc(CFMessagePortRef local, SInt32 msgid, CFDataRef data, void *info)
{
    NSLog(@"this myCallbackFunc");
    return nil;
}

int main(int argc, char *argv[])
{
    CFRunLoopSourceRef source;
    CFRunLoopSourceContext source_context;
    CFRunLoopTimerRef timer;
    CFRunLoopTimerContext timer_context;
    
    bzero(&source_context, sizeof(source_context));
    source_context.perform = _perform;
    source = CFRunLoopSourceCreate(NULL, 0, &source_context);
    CFRunLoopAddSource(CFRunLoopGetCurrent(), source, kCFRunLoopCommonModes);
    
    bzero(&timer_context, sizeof(timer_context));
    timer_context.info = source;
    timer = CFRunLoopTimerCreate(NULL, CFAbsoluteTimeGetCurrent(), 1, 0, 0,
                                 _timer, &timer_context);
    CFRunLoopAddTimer(CFRunLoopGetCurrent(), timer, kCFRunLoopCommonModes);
    
    CFRunLoopRun();
    
    
    return 0;
}


#endif
