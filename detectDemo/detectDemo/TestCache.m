//
//  TestCache.m
//  detectDemo
//
//  Created by renxin on 15/7/6.
//  Copyright (c) 2015年 dhd. All rights reserved.
//

#import "TestCache.h"

@implementation TestCache

+(void) testCache
{
    NSLog(@"This objcet is %p", self);
    NSLog(@"Class is %@, super class is %@", [self class], [self superclass]);
    
    Class currentClass = [self class];
    //
    for (int i = 0; i < 4; i++) {
        NSLog(@"Following the isa pointer %d times gives %p", i, currentClass);
        currentClass = objc_getClass((__bridge void *)currentClass);
    }
    
    NSLog(@"NSObject's class is %p", [NSObject class]);
    NSLog(@"NSObject's meta class is %p", objc_getClass((__bridge void *)[NSObject class]));
}

@end
