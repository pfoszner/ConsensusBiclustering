<?php

$hit_count = @file_get_contents('count-win.txt');
$hit_count++;
@file_put_contents('count-win.txt', $hit_count);

header('Location: http://iiirs.org/IIL/mtba/mtba-win.zip'); // redirect to the real file to be downloaded

?>
